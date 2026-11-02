import numpy as np
import chaospy as cp
import scipy.io as sio
import sklearn.linear_model

def build_design_matrices(order, x, *, with_derivatives=True):
    """Build a PCE basis and design matrices for normalized inputs.

    x: shape (n_samples, n_parameters), with inputs in [0, 1].

    Returns:
        basis: polynomial basis
        Psi: shape (n_samples, n_coefficients)
        Psi_der: shape (n_parameters * n_samples, n_coefficients),
            ordered by parameter, then by sample; None when omitted.

    Derivatives are taken with respect to the normalized inputs.
    """
    x = np.asarray(x)
    n_samples, n_parameters = x.shape

    # Create one independent Uniform(0,1) variable per active parameter
    marginal_distributions = [cp.Uniform(0, 1) for _ in range(n_parameters)]
    # create joint distribution (separate case for n=1)
    if n_parameters == 1:
        joint_distribution = marginal_distributions[0]
    else:
        joint_distribution = cp.J(*marginal_distributions)
    
    basis = cp.generate_expansion(order, joint_distribution, normed=True)
    Psi = cp.call(basis, x.T).T

    # do not compute derivative basis if not required
    if not with_derivatives:
        return basis, Psi, None

    basis_der = cp.gradient(basis)
    Psi_der = cp.call(basis_der, x.T).T
    # Chaospy can omit the parameter axis for one-dimensional inputs.
    if Psi_der.ndim == 2 and n_parameters == 1:
        Psi_der = Psi_der[:, :, None]

    # Parameter 1/all samples, parameter 2/all samples, ...
    Psi_der = Psi_der.transpose(2, 0, 1).reshape(n_parameters * n_samples, Psi.shape[1])

    return basis, Psi, Psi_der

def create_regression_model(regression_type):
    """Create a new, unfitted regression model."""

    model_factories = { "ols": lambda: sklearn.linear_model.LinearRegression(fit_intercept=False),
                        "lasso": lambda: sklearn.linear_model.Lasso(alpha=5e-4,fit_intercept=False,max_iter=30_000,),
                        "omp": lambda: sklearn.linear_model.OrthogonalMatchingPursuit(fit_intercept=False),
                        "lars": lambda: sklearn.linear_model.Lars(fit_intercept=False,n_nonzero_coefs=100),
                        "lassolars": lambda: sklearn.linear_model.LassoLars(fit_intercept=False),
    }

    if regression_type not in model_factories:
        valid_types = [*model_factories.keys()]
        raise ValueError(f"Unknown regression type '{regression_type}'. "f"Choose from {valid_types}.")

    return model_factories[regression_type]()

def calculate_gamma(values, derivatives):
    """Compute a safe gradient scaling factor."""
    derivative_norm = np.linalg.norm(derivatives)
    if derivative_norm <= np.finfo(float).eps:
        raise ValueError("Cannot compute gamma because the derivative norm is zero.")
    return np.linalg.norm(values) / derivative_norm

def normalize_design_columns(design_matrix):
    """Scale every regression column to unit RMS without centering it."""

    design_matrix = np.asarray(design_matrix, dtype=float)

    column_scales = np.linalg.norm(design_matrix, axis=0) / np.sqrt(design_matrix.shape[0])
    if np.any(column_scales <= np.finfo(float).eps):
        zero_columns = np.flatnonzero(column_scales <= np.finfo(float).eps)
        raise ValueError(f"Cannot normalize zero design-matrix columns: {zero_columns.tolist()}")

    return design_matrix / column_scales, column_scales

def fit_independent_models(Psi, Psi_der, y_train_red, y_train_der_red, regression_type, verbosity=False):
    """Fit one standard and one gradient-enhanced model per POD mode."""

    n_coefficients = Psi.shape[1]
    n_modes = y_train_red.shape[1]

    # init coeff matrices
    gepce_coeffs = np.zeros((n_coefficients, n_modes))
    pce_coeffs = np.zeros((n_coefficients, n_modes))
    gamma_values = np.zeros(n_modes)

    # Psi is identical for every POD mode, so its column normalization is
    # computed only once. The GE-PCE matrix is normalized inside the loop
    # because the mode-dependent gamma changes its columns.
    Psi_scaled, pce_column_scales = normalize_design_columns(Psi)

    for mode in range(n_modes):
        if verbosity:
            print(f"Calculating coefficients for POD mode {mode + 1}")
        # compute gamme for every mode
        gamma = calculate_gamma(y_train_red[:, mode],y_train_der_red[:, mode])
        gamma_values[mode] = gamma

        Psi_aug = np.vstack([Psi,gamma * Psi_der,])
        y_aug = np.concatenate([y_train_red[:, mode],gamma * y_train_der_red[:, mode]])

        # column normalization
        Psi_aug_scaled, gepce_column_scales = normalize_design_columns(Psi_aug)

        gepce_regressor = create_regression_model(regression_type)
        pce_regressor = create_regression_model(regression_type)

        gepce_regressor.fit(Psi_aug_scaled, y_aug)
        pce_regressor.fit(Psi_scaled,y_train_red[:, mode])

       # Transform coefficients back to the original polynomial basis.
        gepce_coeffs[:, mode] = gepce_regressor.coef_ / gepce_column_scales
        pce_coeffs[:, mode] = pce_regressor.coef_ / pce_column_scales

    return gepce_coeffs, pce_coeffs, gamma_values

def predict_pod_fields(pce_model, x, U, feature_means, feature_stds):
    """Predict full fields from normalized parameter inputs.

    x: (n_samples, n_parameters)
    U: (n_dofs, n_modes)
    feature_means, feature_stds: (1, n_modes)

    Returns fields shaped (n_dofs, n_samples).
    """
    coefficients_normalized = np.asarray(pce_model(*x.T)).reshape(U.shape[1], x.shape[0])
    coefficients = (coefficients_normalized * feature_stds.T+ feature_means.T)
    return U @ coefficients

def predict_ae_fields(pce_model, x, basis_data):
    """Predict full fields from Autoencoder compressed parameter inputs.

    x: (n_samples, n_parameters)

    Returns fields shaped (n_dofs, n_samples).
    """
    import torch
    n_samples = x.shape[0]
    n_latent = basis_data["feature_means"].shape[1]

    latent_normalized = np.asarray(pce_model(*x.T))
    latent_normalized = latent_normalized.reshape(n_latent, n_samples)
    latent_variables = (latent_normalized* basis_data["feature_stds"].T + basis_data["feature_means"].T).T

    # load AE model
    model = basis_data["model"]
    device = next(model.parameters()).device

    # decode using AE
    latent_tensor = torch.from_numpy(np.ascontiguousarray(latent_variables, dtype=np.float32)).to(device)
    model.eval()
    with torch.no_grad():
        decoded_normalized = model.decode(latent_tensor).cpu().numpy()

    # rescale
    return (basis_data["mean_field"]+ basis_data["field_scale"]* decoded_normalized.T)

def predict_and_save_fields(compression_mode, pce_model, basis_data, output_filepath=None):
    """ Predict fields using the PCE model and save them to a CSV file."""
    # load input basis data (normalized)
    x_train = basis_data["x_train"]
    x_val = basis_data["x_val"]
    x_test = basis_data["x_test"]

    if output_filepath is None:
        raise ValueError("output_filepath must be provided to save the predicted fields.")

    if compression_mode == "POD":
        # predict, denormalize, and reconstruct full fields using inverse POD
        U = basis_data["U"]
        feature_stds = basis_data["feature_stds"]
        feature_means = basis_data["feature_means"]

        # field prediction
        Y_train_pred = Y_train_pred = predict_pod_fields(pce_model, x_train, U, feature_means, feature_stds)
        Y_test_pred = predict_pod_fields(pce_model, x_test, U, feature_means, feature_stds)
        Y_val_pred = predict_pod_fields(pce_model, x_val, U, feature_means, feature_stds)

    elif compression_mode == "None":
        # only predict and denormalize
        feature_stds = basis_data["feature_stds"]
        feature_means = basis_data["feature_means"]
        
        # prediction of normalized fields
        y_train_pred_norm = pce_model(*x_train.T)
        y_test_pred_norm = pce_model(*x_test.T)
        y_val_pred_norm = pce_model(*x_val.T)
        # denormalize fields
        Y_train_pred = y_train_pred_norm * feature_stds + feature_means
        Y_test_pred = y_test_pred_norm * feature_stds + feature_means
        Y_val_pred = y_val_pred_norm * feature_stds + feature_means

    elif compression_mode == "AE":
        Y_train_pred = predict_ae_fields(pce_model, x_train, basis_data)
        Y_test_pred = predict_ae_fields(pce_model, x_test, basis_data)
        Y_val_pred = predict_ae_fields(pce_model, x_val, basis_data)
    else:
        raise ValueError("Choose a suitable compression mode!")

    matlab_data = {
        # Fields use the convention (n_dofs, n_samples)
        "Y_train_pred": Y_train_pred,
        "Y_val_pred": Y_val_pred,
        "Y_test_pred": Y_test_pred,

        # Diagnostic reference data
        "Y_test_true": basis_data["Y_test"],
        "K_mat": basis_data["K_mat"],
    }

    print(f"Saving predicted fields to {output_filepath}...")
    sio.savemat(output_filepath,matlab_data,do_compression=True,long_field_names=True,oned_as="column",)

    return Y_train_pred, Y_val_pred, Y_test_pred