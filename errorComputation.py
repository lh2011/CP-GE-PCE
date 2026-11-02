import numpy as np
from autoencoder import reconstruct_fields
from pod import project_onto_pod_modes
from pce import predict_pod_fields, predict_ae_fields

def compute_errors(compression_mode, pce_model, verbosity, basis_data):
    # unpack basis data required for all cases
    x_train = basis_data["x_train"]
    x_val = basis_data["x_val"]
    x_test = basis_data["x_test"]
    Y_train = basis_data["Y_train"]
    Y_val = basis_data["Y_val"]
    Y_test = basis_data["Y_test"]
    feature_means = basis_data["feature_means"]
    feature_stds = basis_data["feature_stds"]
    K_mat = basis_data["K_mat"]

    if compression_mode == "POD":
        U = basis_data["U"]
        errors_train = calc_error(pce_model, x_train, Y_train, feature_means, feature_stds, U, K_mat)
        errors_val = calc_error(pce_model, x_val, Y_val, feature_means, feature_stds, U, K_mat)
        errors_test = calc_error(pce_model, x_test, Y_test, feature_means, feature_stds, U, K_mat)
    elif compression_mode == "None":
        errors_train = calc_error_no_compression(pce_model, x_train, Y_train, feature_means, feature_stds, K_mat)
        errors_val = calc_error_no_compression(pce_model, x_val, Y_val, feature_means, feature_stds, K_mat)
        errors_test = calc_error_no_compression(pce_model, x_test, Y_test, feature_means, feature_stds, K_mat)
    elif compression_mode == "AE":
        errors_train = calc_error_AE(pce_model, x_train, Y_train, K_mat, basis_data)
        errors_val = calc_error_AE(pce_model, x_val, Y_val, K_mat, basis_data)
        errors_test = calc_error_AE(pce_model, x_test, Y_test, K_mat, basis_data)
    else:
        raise ValueError("Choose a suitable compression mode!")

    errors = {'train': errors_train, 'test': errors_test, 'val': errors_val}

    if verbosity:
        # Print statistics for each error set
        for key in ['train', 'test', 'val']:
            arr = np.array(errors[key])
            print(f"{key.capitalize()} error: mean={arr.mean():.4e}, std={arr.std():.4e}, max={arr.max():.4e}")

    return errors

def calc_k_rel_err(Y_pred,Y,K_mat):
    diff = Y_pred - Y
    diff_K = diff.T @ K_mat
    numerator = np.sum(diff_K * diff.T, axis=1)

    original_K = Y.T @ K_mat  
    denominator = np.sum(original_K * Y.T, axis=1)
    
    errors = np.sqrt(numerator / denominator)
    return errors

def compute_pod_error(verbosity, Y, U, K_mat):
    y_pred_red = project_onto_pod_modes(Y, U, K_mat)
    Y_pred = np.asarray(y_pred_red.T @ U.T).transpose()
    errors = calc_k_rel_err(Y_pred, Y, K_mat)
    if verbosity:
        print(f"\nPOD error: mean={errors.mean():.4e}, std={errors.std():.4e}, max={errors.max():.4e}")
    return errors

def compute_AE_error(verbosity, Y, basis_data):
    K_mat = basis_data["K_mat"]
    model = basis_data["model"]
    mean_field = basis_data["mean_field"]
    field_scale = basis_data["field_scale"]
    # normalize Y
    Y_normed = np.ascontiguousarray(((Y - mean_field) / field_scale).T, dtype=np.float32)
    Y_pred = reconstruct_fields(model, Y_normed, mean_field, field_scale, next(basis_data["model"].parameters()).device)
    errors = calc_k_rel_err(Y_pred, Y, K_mat)
    if verbosity:
        print(f"\nAE reconstruction error: mean={errors.mean():.4e}, std={errors.std():.4e}, max={errors.max():.4e}")
    return errors
    
def calc_error(pce_model, x, Y, feature_means, feature_stds, U, K_mat):
    Y_pred = predict_pod_fields(pce_model, x, U, feature_means, feature_stds)
    return calc_k_rel_err(Y_pred, Y, K_mat)

def calc_error_no_compression(pce_model, x, Y, feature_means, feature_stds, K_mat):
    y_pred_norm = pce_model(*x.T)
    # unnormalize the reduced predictions
    Y_pred = y_pred_norm * feature_stds + feature_means
    errors = calc_k_rel_err(Y_pred, Y, K_mat)
    return errors

def calc_error_AE(pce_model, x, Y, K_mat, basis_data):
    Y_pred = predict_ae_fields(pce_model, x, basis_data)
    errors = calc_k_rel_err(Y_pred, Y, K_mat)
    return errors