import numpy as np
import chaospy as cp
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re

from pce import build_design_matrices, create_regression_model, calculate_gamma, normalize_design_columns 

from pod import get_basis
from dataLoadingHelper import *
from errorComputation import *
from compute_PCE import *

def compute_energy_scores(Y_true, Y_pred, K_mat):
    """
    Y_true, Y_pred: shape (n_dofs, n_samples)
    K_mat: stiffness matrix
    """
    if Y_true.shape != Y_pred.shape:
        raise ValueError("True and predicted fields have different shapes.")

    n_samples = Y_true.shape[1]
    scores = np.empty(n_samples)

    for n in range(n_samples):
        error = Y_true[:, n] - Y_pred[:, n]
        energy_squared = error @ (K_mat @ error)

        # Protect against tiny negative values from roundoff
        scores[n] = np.sqrt(max(energy_squared, 0.0))

    return scores

def conformal_quantile(scores, alpha=0.05):
    scores = np.asarray(scores).reshape(-1)

    if np.any(~np.isfinite(scores)):
        raise ValueError("Scores contain NaN or infinite values.")

    if np.any(scores < 0):
        raise ValueError("Nonconformity scores must be non-negative.")

    N = len(scores)
    k = int(np.ceil((N + 1) * (1.0 - alpha)))

    if k > N:
        return np.inf

    # k-1 because Python uses zero-based indexing
    return np.partition(scores, k - 1)[k - 1]

def conformal_prediction_POD_pce(order, *, num_ds=None, dataset_csv=None, filename=None):
    # parameter definitions
    angle_mode=SAMPLED_ANGLE_MODE
    modes = 32
    partition = 4
    parameter_names = ["alpha","MH","MW","MAG","Theta1"]
    dataset_num_params = 5
    if num_ds is None:
        num_ds = 623
    if filename is None:
        filename = (f"data/basis_p{dataset_num_params}_part{partition}_modes{modes}_aug0.mat")

    # load or compute basis
    basis_data, time_basis = get_basis(filename, partition, parameter_names, modes, use_saved_basis=True, is_augmented=False, dataset_num_params=dataset_num_params, num_ds=num_ds, angle_mode=angle_mode, dataset_csv=dataset_csv)
    # unpack basis data
    x_train = basis_data["x_train"]
    x_val = basis_data["x_val"]
    x_test = basis_data["x_test"]
    Y_train = basis_data["Y_train"]
    Y_val = basis_data["Y_val"]
    Y_test = basis_data["Y_test"]
    y_train_red = basis_data["y_train_red"]
    y_train_der_red = basis_data["y_train_der_red"]
    feature_means = basis_data["feature_means"]
    feature_stds = basis_data["feature_stds"]
    U = basis_data["U"]
    K_mat = basis_data["K_mat"]

    # definitions
    n_modes = y_train_red.shape[1]

    pce_basis, Psi, _ = build_design_matrices(order, x_train, False)
    n_coefficients = Psi.shape[1]

    pce_coeffs = np.zeros((n_coefficients, n_modes))

    # obtain pce coeffs for every mode
    Psi_scaled, pce_column_scales = normalize_design_columns(Psi)
    for mode in range(n_modes): 
        # create and fit regression model
        pce_regressor = create_regression_model("ols")
        pce_regressor.fit(Psi_scaled,y_train_red[:, mode])

        # Transform coefficients back to the original polynomial basis.
        pce_coeffs[:, mode] = pce_regressor.coef_ / pce_column_scales

    # build model from coeffs
    pce_model = cp.sum(pce_basis[:, None] * pce_coeffs,axis=0)

    # get all training predictions and errors
    y_train_pred_red = pce_model(*x_train.T).T
    E_train = y_train_red - y_train_pred_red

    # Leverage values
    # this is a nicer way to compute h_mm
    Q, _ = np.linalg.qr(Psi, mode="reduced")
    leverage = np.sum(Q**2, axis=1)

    # LOO predictions of normalized POD coefficients
    y_train_pred_red_loo = y_train_red - E_train / (1.0 - leverage[:, None])

    # Reverse normalization
    y_train_pred_loo = y_train_pred_red_loo * feature_stds + feature_means

    # Reconstruct complete LOO fields
    Y_train_pred_loo = U @ y_train_pred_loo.T

    scores_train = compute_energy_scores(Y_train, Y_train_pred_loo, K_mat)

    return scores_train, pce_model, basis_data

def conformal_prediction_POD_gepce(order, *, num_ds, dataset_csv=None, filename=None):
    # parameter definitions
    angle_mode=SAMPLED_ANGLE_MODE
    modes = 32
    partition = 4
    parameter_names = ["alpha","MH","MW","MAG","Theta1"]
    dataset_num_params = 5
    if num_ds is None:
        num_ds = 623
    if filename is None:
        filename = (f"data/basis_p{dataset_num_params}_part{partition}_modes{modes}_aug0.mat")

    # load or compute basis
    basis_data, time_basis = get_basis(filename, partition, parameter_names, modes, use_saved_basis=True, is_augmented=False, dataset_num_params=dataset_num_params, num_ds=num_ds, angle_mode=angle_mode, dataset_csv=dataset_csv)
    # unpack basis data
    x_train = basis_data["x_train"]
    x_val = basis_data["x_val"]
    x_test = basis_data["x_test"]
    Y_train = basis_data["Y_train"]
    Y_val = basis_data["Y_val"]
    Y_test = basis_data["Y_test"]
    y_train_red = basis_data["y_train_red"]
    y_train_der_red = basis_data["y_train_der_red"]
    feature_means = basis_data["feature_means"]
    feature_stds = basis_data["feature_stds"]
    U = basis_data["U"]
    K_mat = basis_data["K_mat"]

    # definitions
    n_parameters = len(parameter_names)
    n_train = x_train.shape[0]
    n_test = x_test.shape[0]
    n_modes = y_train_red.shape[1]

    pce_basis, Psi, Psi_der = build_design_matrices(order, x_train)
    n_coefficients = Psi.shape[1]

    # this is stupid but I need to obtain the blockwise formulation and I dont want to optimize right now
    Psi_der_blocks = Psi_der.reshape(n_parameters, n_train, n_coefficients).transpose(1, 0, 2)
    # shape: (n_train, n_parameters, n_coefficients)
    y_der_blocks = y_train_der_red.reshape(n_parameters, n_train, n_modes).transpose(1, 0, 2)

    # values required for looping
    n_coefficients = Psi.shape[1]
    gamma_values = np.zeros(n_modes)
    gepce_coeffs = np.zeros((n_coefficients, n_modes))
    # for LOO residual calculation
    loo_residuals_red = np.empty((n_train, n_modes))
    block_size = n_parameters + 1

    for mode in range(n_modes):
        # compute gamme for every mode
        gamma = calculate_gamma(y_train_red[:, mode],y_train_der_red[:, mode])
        gamma_values[mode] = gamma

        Psi_aug_blocks = np.concatenate([Psi[:, None, :], gamma * Psi_der_blocks],axis=1)
        # Observation blocks for this mode:
        # shape (n_train, 1+d)
        y_aug_blocks = np.concatenate([y_train_red[:, mode, None], gamma * y_der_blocks[:, :, mode]],axis=1)

        # Flatten the sample-major blocks for regression
        Psi_aug = Psi_aug_blocks.reshape(n_train * (n_parameters + 1),n_coefficients,)
        y_aug = y_aug_blocks.reshape(n_train * (n_parameters + 1))

        # column normalization
        Psi_aug_scaled, gepce_column_scales = normalize_design_columns(Psi_aug)
        # regression using OLS
        gepce_regressor = create_regression_model("ols")
        gepce_regressor.fit(Psi_aug_scaled, y_aug)
        # Transform coefficients back to the original polynomial basis.
        gepce_coeffs[:, mode] = gepce_regressor.coef_ / gepce_column_scales

        # -------------------------------------- LOO residual computation --------------------------------------
        # Ordinary residuals of the augmented regression
        y_aug_pred = gepce_regressor.predict(Psi_aug_scaled)
        residual_aug = y_aug - y_aug_pred

        # Each row now corresponds to one physical sample block:
        # [value residual, derivative residual 1, ..., derivative residual d]
        residual_blocks = residual_aug.reshape(n_train, block_size)
        identity_block = np.eye(block_size)

        # Thin QR decomposition of the same design matrix used for fitting
        Q, R = np.linalg.qr(Psi_aug_scaled, mode="reduced")
        if np.linalg.matrix_rank(R) < n_coefficients:
            raise np.linalg.LinAlgError(f"GE-PCE design matrix for mode {mode} is rank deficient.")
        # Rows corresponding to each physical sample
        Q_blocks = Q.reshape(n_train, block_size, Q.shape[1])

        # compute residual for every training sample
        for n in range(n_train):
            Q_n = Q_blocks[n]

            # Block leverage matrix, shape (d+1, d+1)
            H_nn = Q_n @ Q_n.T

            # Complete block LOO residual
            loo_residual_block = np.linalg.solve(identity_block - H_nn,residual_blocks[n],)

            # The first row is the function-value residual
            # this corresponds to the error selector e_0
            loo_residuals_red[n, mode] = loo_residual_block[0]

    # build model from coeffs
    gepce_model = cp.sum(pce_basis[:, None] * gepce_coeffs,axis=0)

    # Normalized reduced LOO predictions
    y_train_pred_red_loo = y_train_red - loo_residuals_red

    # Reverse POD-coefficient normalization
    y_train_pred_loo = y_train_pred_red_loo * feature_stds+ feature_means
    
    # Reconstruct physical fields
    # Shape: (n_dofs, n_train)
    Y_train_pred_loo = U @ y_train_pred_loo.T

    # Compute the energy-norm LOO scores
    scores_train = compute_energy_scores(Y_train, Y_train_pred_loo, K_mat)

    return scores_train, gepce_model, basis_data

def pce_test(order, nominal_coverages, *, num_ds = None, dataset_csv=None, filename=None):
    scores_train, pce_model, basis_data = conformal_prediction_POD_pce(order, num_ds=num_ds, dataset_csv=dataset_csv, filename=filename)

    # unpack basis data
    x_test = basis_data["x_test"]
    Y_test = basis_data["Y_test"]
    feature_means = basis_data["feature_means"]
    feature_stds = basis_data["feature_stds"]
    U = basis_data["U"]
    K_mat = basis_data["K_mat"]

    # predict test data
    y_test_pred_red = pce_model(*x_test.T).T
    y_test_pred = y_test_pred_red * feature_stds + feature_means
    Y_test_pred = U @ y_test_pred.T

    # compute test scores
    scores_test = compute_energy_scores(Y_test, Y_test_pred, K_mat)
    coverages = []
    radii = []

    for nominal_coverage in nominal_coverages:
        alpha_now = 1.0 - nominal_coverage
        q_now = conformal_quantile(scores_train, alpha_now)
        coverage_now = np.mean(scores_test <= q_now)
        coverages.append(coverage_now)
        radii.append(q_now)
        print(
            f"Nominal: {100*nominal_coverage:5.1f}% | "
            f"Empirical: {100*coverage_now:5.1f}% | "
            f"Radius: {q_now:.6e}"
        )

    pce_errors = compute_errors("POD", pce_model, verbosity=True, basis_data=basis_data)

    return coverages, radii

def gepce_test(order, nominal_coverages, *, num_ds=None, dataset_csv=None, filename=None):
    scores_train, gepce_model, basis_data = conformal_prediction_POD_gepce(order, num_ds=num_ds, dataset_csv=dataset_csv, filename=filename)
    # unpack basis data
    x_test = basis_data["x_test"]
    Y_test = basis_data["Y_test"]
    feature_means = basis_data["feature_means"]
    feature_stds = basis_data["feature_stds"]
    U = basis_data["U"]
    K_mat = basis_data["K_mat"]
    
    # predict test data
    y_test_pred_red = gepce_model(*x_test.T).T
    y_test_pred = y_test_pred_red * feature_stds + feature_means
    Y_test_pred = U @ y_test_pred.T
    
    # compute test scores
    scores_test = compute_energy_scores(Y_test, Y_test_pred, K_mat)
    coverages = []
    radii = []
    
    for i,nominal_coverage in enumerate(nominal_coverages):
        alpha_now = 1.0 - nominal_coverage
        q_now = conformal_quantile(scores_train, alpha_now)
        coverage_now = np.mean(scores_test <= q_now)
        coverages.append(coverage_now)
        radii.append(q_now)
    
        print(
            f"Nominal: {100*nominal_coverage:5.1f}% | "
            f"Empirical: {100*coverage_now:5.1f}% | "
            f"Radius: {q_now:.6e}"
        )

    gepce_errors = compute_errors("POD", gepce_model, verbosity=True, basis_data=basis_data)

    return coverages, radii

def run_test_for_multiple_datasets(order=3,nominal_coverages=[0.80, 0.90, 0.95, 0.99],data_dir="data/cp_eval_data"):
    data_dir = Path(data_dir)
    pattern = re.compile(r"dataset_p_5_samp_(?P<num_ds>\d+)_seed_(?P<seed>\d+)\.csv")

    # extract all file names, including seed and number of samples
    datasets = []
    for csv_file in data_dir.glob("dataset_p_5_samp_*_seed_*.csv"):
        match = pattern.fullmatch(csv_file.name)
        if match:
            datasets.append((
                int(match.group("seed")),
                int(match.group("num_ds")),
                csv_file,
            ))
    datasets.sort()

    results = []
    output_file = data_dir / f"conformal_results_order_{order}.csv"

    for seed, num_ds, csv_file in datasets:
        print(f"Computing data for seed {seed}...")

        basis_file = data_dir / f"basis_p5_part4_modes32_aug0_samp_{num_ds}_seed_{seed}.mat"

        # compute pce values
        print(csv_file)
        pce_coverages, pce_radii = pce_test(order, nominal_coverages, num_ds=num_ds, dataset_csv=csv_file, filename=basis_file)
        # compute gepce values
        gepce_coverages, gepce_radii = gepce_test(order, nominal_coverages, num_ds=num_ds, dataset_csv=csv_file, filename=basis_file)

        # this saves all the data for one seed
        for method, coverages, radii in (("PCE", pce_coverages, pce_radii), ("GE-PCE", gepce_coverages, gepce_radii)):
                # this writes one row
                for nominal, empirical, radius in zip(nominal_coverages, coverages, radii, strict=True):
                    results.append({
                        "seed": seed,
                        "num_ds": num_ds,
                        "method": method,
                        "nominal_coverage": nominal,
                        "empirical_coverage": empirical,
                        "radius": radius,
                    })

        # Save after each completed seed so earlier results survive a later error.
        pd.DataFrame(results).to_csv(output_file, index=False)

    return pd.DataFrame(results)

def visualize_tests(nominal_coverages, order):

    df = pd.read_csv(f"data/cp_eval_data/conformal_results_order_{order}.csv")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    colors = {"PCE": "tab:blue", "GE-PCE": "tab:orange"}

    # 1. Each seed's calibration curve, with the mean emphasized
    for method, color in colors.items():
        subset = df[df["method"] == method]

        for _, seed_data in subset.groupby("seed"):
            seed_data = seed_data.sort_values("nominal_coverage")
            ax1.plot(
                seed_data["nominal_coverage"],
                seed_data["empirical_coverage"],
                color=color, alpha=0.20, linewidth=1,
            )

        mean = subset.groupby("nominal_coverage")["empirical_coverage"].mean()
        ax1.plot(mean.index, mean.values, "o-",
                color=color, linewidth=2.5, label=f"{method} mean")

    ax1.plot([0, 1], [0, 1], "k--", label="Nominal")
    ax1.set(xlabel="Nominal coverage", ylabel="Empirical coverage",
            xlim=(min(nominal_coverages)-0.02, 1.0), ylim=(min(nominal_coverages)-0.02, 1.01))
    ax1.legend()
    ax1.set_title(f"Empirical vs. Nominal Coverage (Order {order})")
    ax1.grid(alpha=0.3)

    # 2. Paired results for each seed at 95%
    at95 = df[np.isclose(df["nominal_coverage"], 0.95)]
    wide = at95.pivot(index="seed", columns="method",
                    values="empirical_coverage").sort_index()

    x = np.arange(len(wide))
    for i, (pce, gepce) in enumerate(zip(wide["PCE"], wide["GE-PCE"])):
        ax2.plot(
            [i - 0.12, i + 0.12],
            [pce, gepce],
            color="0.7",
            linewidth=1,
        )

    for method, offset in (("PCE", -0.12), ("GE-PCE", 0.12)):
        ax2.scatter(x + offset, wide[method],
                    color=colors[method], label=method, zorder=3)

    ax2.axhline(0.95, color="black", linestyle="--", label="Nominal 95%")
    ax2.set(xlabel="Seed", ylabel="Empirical coverage at 95%")
    ax2.set_xticks(x, wide.index)
    ax2.grid(axis="y", alpha=0.3)
    ax2.set_title(f"Coverage at 0.95 Target Across Seeds (Order {order})")
    ax2.legend()

    fig.tight_layout()
    fig.savefig(f"data/cp_eval_data/coverage_across_seeds_order_{order}.png",
                dpi=200, bbox_inches="tight")

    radii = df.pivot(
        index=["seed", "nominal_coverage"],
        columns="method",
        values="radius",
    ).reset_index()

    radii["GE/PCE"] = radii["GE-PCE"] / radii["PCE"]

    fig_radii, ax_radii = plt.subplots(figsize=(7, 4.5))

    for seed, group in radii.groupby("seed"):
        group = group.sort_values("nominal_coverage")
        ax_radii.plot(group["nominal_coverage"], group["GE/PCE"],
                "o-", color="tab:blue", alpha=0.25)

    median = radii.groupby("nominal_coverage")["GE/PCE"].median()
    ax_radii.plot(median.index, median.values, "o-",
            color="tab:blue", linewidth=3, label="Median across seeds")

    ax_radii.axhline(1, color="black", linestyle="--", label="Equal radii")
    ax_radii.set(xlabel="Nominal coverage", ylabel="GE-PCE radius / PCE radius")
    ax_radii.grid(alpha=0.3)
    ax_radii.set_title(f"Radii for 95% Nominal Coverage (Order {order})")
    ax_radii.legend()

    fig_radii.savefig(f"data/cp_eval_data/radii_for_095_order_{order}.png",
        dpi=200, bbox_inches="tight")

    plt.show()

def main():
    order = 3
    nominal_coverages = [0.80, 0.90, 0.95, 0.99]

    results = run_test_for_multiple_datasets(order=order,nominal_coverages=nominal_coverages,data_dir="data/cp_eval_data")
    visualize_tests(nominal_coverages, order)
    #print("Standard POD-PCE:")
    #pce_test(order, nominal_coverages)
    #print("Gradient-Enhanced POD-PCE:")
    #gepce_test(order, nominal_coverages)

if __name__ == "__main__":
    main()
