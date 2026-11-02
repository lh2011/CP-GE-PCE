import numpy as np
import chaospy as cp
import pandas as pd
import matplotlib.pyplot as plt
import time
from pathlib import Path

from parameters import ALL_PARAMETERS
from pce import build_design_matrices, fit_independent_models
from pod import get_basis, perform_pod
from autoencoder import load_autoencoder_checkpoint, reconstruct_fields, relative_energy_errors, get_basis_AE
from dataLoadingHelper import *
from errorComputation import *

VALID_REGRESSION_TYPES = {
        "ols",
        "lasso",
        "omp",
        "lars",
        "lassolars",
    }

def gepce(modes,order,partition,parameter_names,regression_type,filename="",is_augmented=False,use_saved_basis=False,verbosity=False, angle_mode=LEGACY_ANGLE_MODE, *, num_ds=163, dataset_num_params=None):
    # error checking for input arguments
    parameter_names = list(parameter_names)
    if not parameter_names:
        raise ValueError("At least one parameter must be selected.")
    if len(parameter_names) != len(set(parameter_names)):
        raise ValueError("parameter_names contains duplicates.")
    if regression_type not in VALID_REGRESSION_TYPES:
        raise ValueError(f"Unknown regression type '{regression_type}'. "f"Choose from {sorted(VALID_REGRESSION_TYPES)}.")
    
    # load or compute basis
    basis_data, time_basis = get_basis(filename, partition, parameter_names, modes, use_saved_basis, is_augmented, dataset_num_params=dataset_num_params, num_ds=num_ds, angle_mode=angle_mode)
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

    start = time.time()
    # check for correct shapes
    n_parameters = len(parameter_names)
    n_samples = x_train.shape[0]
    n_modes = y_train_red.shape[1]

    if x_train.shape[1] != n_parameters:
        raise ValueError(f"x_train contains {x_train.shape[1]} parameters, but {n_parameters} parameter names were provided.")
    if n_modes != modes:
        raise ValueError(f"The reduced data contains {n_modes} modes, but modes={modes} was requested.")
    expected_derivative_rows = n_parameters * n_samples
    if y_train_der_red.shape != (expected_derivative_rows,n_modes,):
        raise ValueError(f"Unexpected derivative-data shape. Expected {(expected_derivative_rows, n_modes)}, but received {y_train_der_red.shape}.")

    pce_basis, Psi, Psi_der = build_design_matrices(order, x_train)

    # fit indepenedent model for each mode using 
    gepce_coeffs, pce_coeffs, gamma_values = fit_independent_models(Psi=Psi,Psi_der=Psi_der,y_train_red=y_train_red,y_train_der_red=y_train_der_red,regression_type=regression_type, verbosity=verbosity)

    # Construct vector-valued polynomial models
    gepce_model = cp.sum(pce_basis[:, None] * gepce_coeffs,axis=0)
    pce_model = cp.sum(pce_basis[:, None] * pce_coeffs,axis=0)

    end = time.time()
    time_taken = end - start + time_basis

    if verbosity:
        print(f"\nStandard PCE: order={order}, modes={modes}, "f"parameters={parameter_names}, solver={regression_type}")
    pce_errors = compute_errors("POD", pce_model, verbosity, basis_data)

    if verbosity:
        print(f"\nGE-PCE: order={order}, modes={modes}, "f"parameters={parameter_names}, solver={regression_type}")
    gepce_errors = compute_errors("POD", gepce_model, verbosity, basis_data)
    
    if gamma_values is not None and verbosity:
            print("Gradient scaling factors: "f"min={gamma_values.min():.4e}, "f"mean={gamma_values.mean():.4e}, "f"max={gamma_values.max():.4e}")

    errors_pod = compute_pod_error(verbosity, Y_test, U, K_mat)
    sparsity_tolerance = 1e-8
    active_gepce_coeffs = np.count_nonzero(np.abs(gepce_coeffs) >= sparsity_tolerance)
    active_pce_coeffs = np.count_nonzero(np.abs(pce_coeffs) >= sparsity_tolerance)

    return gepce_model, pce_model, gepce_errors, pce_errors, errors_pod, active_gepce_coeffs, active_pce_coeffs, basis_data, time_taken

def gepce_AE(modes,order,partition,parameter_names,regression_type,verbosity=False, angle_mode=LEGACY_ANGLE_MODE, *, checkpoint_path, num_ds=163, dataset_num_params=None):
    # error checking for input arguments
    parameter_names = list(parameter_names)
    if not parameter_names:
        raise ValueError("At least one parameter must be selected.")
    if len(parameter_names) != len(set(parameter_names)):
        raise ValueError("parameter_names contains duplicates.")
    if regression_type not in VALID_REGRESSION_TYPES:
        raise ValueError(f"Unknown regression type '{regression_type}'. "f"Choose from {sorted(VALID_REGRESSION_TYPES)}.")
    
    # load or compute basis
    basis_data, time_basis = get_basis_AE(partition, parameter_names, modes, dataset_num_params=dataset_num_params, checkpoint_path=checkpoint_path, num_ds=num_ds, angle_mode=angle_mode)

    x_train = basis_data["x_train"]
    y_train_red = basis_data["y_train_red"]
    y_train_der_red = basis_data["y_train_der_red"]

    start = time.time()
    # check for correct shapes
    n_parameters = len(parameter_names)
    n_samples = x_train.shape[0]
    n_modes = y_train_red.shape[1]

    if x_train.shape[1] != n_parameters:
        raise ValueError(f"x_train contains {x_train.shape[1]} parameters, but {n_parameters} parameter names were provided.")
    if n_modes != modes:
        raise ValueError(f"The reduced data contains {n_modes} modes, but modes={modes} was requested.")
    expected_derivative_rows = n_parameters * n_samples
    if y_train_der_red.shape != (expected_derivative_rows,n_modes,):
        raise ValueError(f"Unexpected derivative-data shape. Expected {(expected_derivative_rows, n_modes)}, but received {y_train_der_red.shape}.")

    pce_basis, Psi, Psi_der = build_design_matrices(order, x_train)

    # fit indepenedent model for each mode using 
    gepce_coeffs, pce_coeffs, gamma_values = fit_independent_models(Psi=Psi,Psi_der=Psi_der,y_train_red=y_train_red,y_train_der_red=y_train_der_red,regression_type=regression_type, verbosity=verbosity)

    # Construct vector-valued polynomial models
    gepce_model = cp.sum(pce_basis[:, None] * gepce_coeffs,axis=0)
    pce_model = cp.sum(pce_basis[:, None] * pce_coeffs,axis=0)

    end = time.time()
    time_taken = end - start + time_basis

    if verbosity:
        print(f"\nStandard PCE: order={order}, modes={modes}, "f"parameters={parameter_names}, solver={regression_type}")
    pce_errors = compute_errors("AE", pce_model, verbosity, basis_data)

    if verbosity:
        print(f"\nGE-PCE: order={order}, modes={modes}, "f"parameters={parameter_names}, solver={regression_type}")
    gepce_errors = compute_errors("AE", gepce_model, verbosity, basis_data)
    if gamma_values is not None and verbosity:
            print("Gradient scaling factors: "f"min={gamma_values.min():.4e}, "f"mean={gamma_values.mean():.4e}, "f"max={gamma_values.max():.4e}")

    errors_AE = compute_AE_error(verbosity, basis_data["Y_test"], basis_data)
    sparsity_tolerance = 1e-8
    active_gepce_coeffs = np.count_nonzero(np.abs(gepce_coeffs) >= sparsity_tolerance)
    active_pce_coeffs = np.count_nonzero(np.abs(pce_coeffs) >= sparsity_tolerance)

    return gepce_model, pce_model, gepce_errors, pce_errors, errors_AE, active_gepce_coeffs, active_pce_coeffs, basis_data, time_taken

def compare_gepce(orders, modes, parameter_names, partitions, num_ds, regression_type, angle_mode):
    dataset_num_params = len(parameter_names)
    is_augmented = False

    csv_rows = []

    for j, order in enumerate(orders):

        err = np.zeros(len(partitions))
        err_normal = np.zeros(len(partitions))
        err_train = np.zeros(len(partitions))
        err_train_normal = np.zeros(len(partitions))
        err_val = np.zeros(len(partitions))
        err_val_normal = np.zeros(len(partitions))
        samples = np.zeros(len(partitions))
        err_pod = np.zeros(len(partitions))

        for i, partition in enumerate(partitions):
            print(f"Computing for order {order} and partition {partition}")
            filename = (f"data/basis_p{dataset_num_params}_part{partition}_modes{modes}_aug{int(is_augmented)}.mat")
            (pce_model,pce_model_normal,errors,errors_normal,errors_pod_test,active_coeffs,active_coeffs_normal, basis_data, time_taken) = gepce(modes,order,partition,parameter_names,regression_type=regression_type,filename=filename,is_augmented=False,use_saved_basis=True,verbosity=False,angle_mode=angle_mode,num_ds=num_ds)

            err[i] = errors["test"].mean()
            err_normal[i] = errors_normal["test"].mean()
            err_train[i] = errors["train"].mean()
            err_train_normal[i] = errors_normal["train"].mean()
            err_val[i] = errors["val"].mean()
            err_val_normal[i] = errors_normal["val"].mean()
            samples[i] = errors["train"].shape[0]
            err_pod[i] = errors_pod_test.mean()

            # Save one row for the CSV
            csv_rows.append({
                "order": order,
                "partition": partition,
                "samples": samples[i],
                "gepce_error": err[i],
                "gfpce_error": err_normal[i],
                "pod_error": err_pod[i],
                "active_ge": active_coeffs,
                "active_gf": active_coeffs_normal,
                "gepce_error_train": err_train[i],
                "gfpce_error_train": err_train_normal[i],
                "gepce_error_val": err_val[i],
                "gfpce_error_val": err_val_normal[i],
            })

    # Write CSV
    df = pd.DataFrame(csv_rows)
    for order in df["order"].unique():
        df[df["order"] == order].to_csv(f"compare_gepce_order_{order}_{regression_type}_modes_{modes}.csv",index=False)


def compare_num_modes_POD_AE():
    """Study POD vs AE."""
    partition = 4
    modes_list = np.arange(4, 61, 4)
    parameter_names = ["alpha", "MH", "MW", "MAG", "Theta1"]
    dataset_num_params = 5
    num_ds = 623

    x_train,x_val,x_test,data_train,ders_train,data_val,data_test,K_mat,ids_train,ids_val,ids_test = load_sampled_alpha_data(num_ds=num_ds, partition=partition, parameter_names=parameter_names, dataset_num_params=dataset_num_params)

    def column_K_norm(values):
        """Return the K-norm of every column in values."""
        squared_norms = np.sum(values * (K_mat @ values), axis=0)
        return np.sqrt(np.maximum(squared_norms, 0.0))

    def reconstruction_residual(values, basis):
        """Return values minus their K-orthogonal POD projection."""
        reduced_values = basis.T @ (K_mat @ values)
        return values - basis @ reduced_values

    # The POD modes are nested, so we compute them only once. Using the first m columns gives the same basis as recomputing POD separately for m modes.
    U_all, _ = perform_pod(data_train, K_mat, int(modes_list.max()))

    train_norms = column_K_norm(data_train)
    test_norms = column_K_norm(data_test)
    if np.any(train_norms == 0.0) or np.any(test_norms == 0.0):
        raise ValueError("A field snapshot has zero K-norm.")

    # How strongly do the test fields vary around the training mean field?
    mean_field = np.mean(data_train, axis=1, keepdims=True)
    fluctuation_norms = column_K_norm(data_test - mean_field)
    relative_variation = fluctuation_norms / test_norms

    print(f"Angle-dependent variation relative to the full field: mean={100*np.mean(relative_variation):.4f}%, max={100*np.max(relative_variation):.4f}%")

    mean_train_error = np.zeros(len(modes_list))
    mean_test_error = np.zeros(len(modes_list))
    max_test_error = np.zeros(len(modes_list))

    mean_train_error_ae = np.zeros(len(modes_list))
    mean_test_error_ae = np.zeros(len(modes_list))
    max_test_error_ae = np.zeros(len(modes_list))

    for j, modes in enumerate(modes_list):
        basis = U_all[:, :modes]

        train_residual = reconstruction_residual(data_train, basis)
        test_residual = reconstruction_residual(data_test, basis)

        train_error = column_K_norm(train_residual) / train_norms
        test_error = column_K_norm(test_residual) / test_norms

        mean_train_error[j] = np.mean(train_error)
        mean_test_error[j] = np.mean(test_error)
        max_test_error[j] = np.max(test_error)

        print(
            f"{modes:2d} modes: "
            f"train mean={100*mean_train_error[j]:.4f}%, "
            f"test mean={100*mean_test_error[j]:.4f}%, "
            f"test max={100*max_test_error[j]:.4f}%, "
        )

        #load autoencoder
        loaded_ae = load_autoencoder_checkpoint(f"models/latent_sweep/ae_p5_r{modes}_h256_128_silu_seed42.pt")

        model = loaded_ae["model"]
        mean_field = loaded_ae["mean_field"]
        field_scale = loaded_ae["field_scale"]
        device = loaded_ae["device"]

        network_train = np.ascontiguousarray(((data_train - mean_field) / field_scale).T, dtype=np.float32)
        network_test = np.ascontiguousarray(((data_test - mean_field) / field_scale).T, dtype=np.float32)

        # reconstruct training and testing fields and compute relative energy errors
        reconstruction_train = reconstruct_fields(model=model, normalized_fields=network_train, mean_field=mean_field, field_scale=field_scale, device=device)
        reconstruction_test = reconstruct_fields(model=model, normalized_fields=network_test, mean_field=mean_field, field_scale=field_scale, device=device)
        
        errors_train_ae = relative_energy_errors(reference_fields=data_train, reconstructed_fields=reconstruction_train, K_mat=K_mat)
        errors_test_ae = relative_energy_errors(reference_fields=data_test, reconstructed_fields=reconstruction_test, K_mat=K_mat)

        mean_train_error_ae[j] = np.mean(errors_train_ae)
        mean_test_error_ae[j] = np.mean(errors_test_ae)
        max_test_error_ae[j] = np.max(errors_test_ae)

        print(
            f"{modes:2d} latent variables: "
            f"train mean={100*mean_train_error_ae[j]:.4f}%, "
            f"test mean={100*mean_test_error_ae[j]:.4f}%, "
            f"test max={100*max_test_error_ae[j]:.4f}%, "
        )

    n_train = len(ids_train)
    n_test = len(ids_test)

    pd.DataFrame({
        "modes": modes_list,
        "pod_train_mean": mean_train_error,
        "pod_test_mean": mean_test_error,
        "pod_test_max": max_test_error,
        "ae_train_mean": mean_train_error_ae,
        "ae_test_mean": mean_test_error_ae,
        "ae_test_max": max_test_error_ae,
    }).to_csv(f"compression_errors_p{dataset_num_params}.csv", index=False)

    fig, axes = plt.subplots(1, 1, figsize=(11, 4))

    axes.semilogy(modes_list, 100 * mean_train_error, "+-", label="POD training mean")
    axes.semilogy(modes_list, 100 * mean_test_error, "*-", label="POD testing mean")
    axes.semilogy(modes_list, 100 * max_test_error, "o-", label="POD testing maximum")

    axes.semilogy(modes_list, 100 * mean_train_error_ae, "+-", label="AE training mean")
    axes.semilogy(modes_list, 100 * mean_test_error_ae, "*-", label="AE testing mean")
    axes.semilogy(modes_list, 100 * max_test_error_ae, "o-", label="AE testing maximum")

    axes.set_xlabel("number of POD modes/latent variables")
    axes.set_ylabel("relative $K$-norm error (%)")
    axes.grid(True, which="both")
    axes.legend()
    axes.set_title(f"Parameter names: {', '.join(parameter_names)}, training samples {n_train}, testing samples {n_test}")
    plt.show()

    return modes_list,mean_train_error,mean_test_error,max_test_error

def gepce_no_POD(order,regression_type,num_ds, partition, parameter_names, dataset_num_params=None, *, angle_mode=SAMPLED_ANGLE_MODE, verbosity=False):
    # load data from files
    x_train, x_val, x_test, data_train, ders_train, data_val, data_test, K_mat, ids_train, ids_val, ids_test  = load_data(
        num_ds=num_ds, partition=partition, parameter_names=parameter_names, dataset_num_params=dataset_num_params, angle_mode=angle_mode)

    start = time.time()

    n_parameters = len(parameter_names)
    n_samples = x_train.shape[0]
    n_dofs = data_train.shape[0]

    if data_train.shape[1] != n_samples:
        raise ValueError(f"Inconsistent sample counts: x_train={x_train.shape}, data_train={data_train.shape}")

    # One mean and standard deviation per DOF, calculated across samples.
    feature_means = np.mean(data_train, axis=1, keepdims=True)  # (dofs, 1)
    feature_stds = np.std(data_train, axis=1, keepdims=True)    # (dofs, 1)

    # Function-value RHS:
    # (DOFs, samples) -> (samples, DOFs)
    y_train = ((data_train - feature_means) / feature_stds).T
    print("y_train:", y_train.shape)

    # Raw derivative arrangement:
    # (DOFs, parameters, samples)
    #
    # Desired derivative RHS arrangement before flattening:
    # (parameters, samples, DOFs)
    derivatives_parameter_first = ders_train.transpose(1, 2, 0)
    # Match the parameter-major ordering used below for Psi_der:
    derivative_rhs = derivatives_parameter_first.reshape(n_parameters * n_samples, n_dofs)
    y_train_der = derivative_rhs / feature_stds.T
    print("y_train_der:", y_train_der.shape)

    pce_basis, Psi, Psi_der = build_design_matrices(order, x_train)

    gepce_coeffs, pce_coeffs, gamma_values = fit_independent_models(Psi=Psi,Psi_der=Psi_der,y_train_red=y_train,y_train_der_red=y_train_der,regression_type=regression_type, verbosity=False)

    # Construct vector-valued polynomial models
    gepce_model = cp.sum(pce_basis[:, None] * gepce_coeffs,axis=0)
    pce_model = cp.sum(pce_basis[:, None] * pce_coeffs,axis=0)

    end = time.time()
    time_taken = end - start

    # annoying but required for the error compuation in the end
    basis_data = {
        "x_train": x_train,
        "x_val": x_val,
        "x_test": x_test,
        "Y_train": data_train,
        "Y_val": data_val,
        "Y_test": data_test,
        "feature_means": feature_means,
        "feature_stds": feature_stds,
        "K_mat": K_mat,
    }

    if verbosity:
        print(f"\nStandard PCE: order={order}, "f"parameters={parameter_names}, solver={regression_type}")
    pce_errors = compute_errors("None", pce_model, verbosity, basis_data)

    if verbosity:
        print(f"\nGE-PCE: order={order}, "f"parameters={parameter_names}, solver={regression_type}")
    gepce_errors = compute_errors("None", gepce_model, verbosity, basis_data)

    return pce_errors, gepce_errors, gamma_values, basis_data, time_taken

def compare_POD_AE(order=4):
    modes_list = np.arange(4, 61, 4)
    partition = 3
    parameter_names = ["alpha", "MH","MW","MAG","Theta1",]
    dataset_num_params = 5
    regression_type = "ols"
    angle_mode=SAMPLED_ANGLE_MODE
    num_ds = 623

    # compute no POD once as baseline and then compute with POD for different numbers of modes
    pce_errors, gepce_errors, gamma_values, basis_data, time_taken_no_POD = gepce_no_POD(order=order,regression_type=regression_type,num_ds=num_ds, partition=partition, parameter_names=parameter_names, dataset_num_params=dataset_num_params)

    n_train = len(gepce_errors["train"])
    n_test = len(gepce_errors["test"])

    err_no_POD = gepce_errors["test"].mean()
    err_normal_no_POD = pce_errors["test"].mean()

    err_no_POD_max = gepce_errors["test"].max()
    err_normal_no_POD_max = pce_errors["test"].max()

    err = np.zeros(len(modes_list))
    err_normal = np.zeros(len(modes_list))
    err_train = np.zeros(len(modes_list))
    err_train_normal = np.zeros(len(modes_list))
    err_max = np.zeros(len(modes_list))
    err_normal_max = np.zeros(len(modes_list))
    err_train_max = np.zeros(len(modes_list))
    err_train_normal_max = np.zeros(len(modes_list))

    err_pod = np.zeros(len(modes_list))
    err_pod_max = np.zeros(len(modes_list))

    times_POD = np.zeros(len(modes_list))
    times_AE = np.zeros(len(modes_list))

    err_ae = np.zeros(len(modes_list))
    err_ae_max = np.zeros(len(modes_list))
    err_train_ae = np.zeros(len(modes_list))
    err_train_ae_max = np.zeros(len(modes_list))

    err_ae_pce = np.zeros(len(modes_list))
    err_ae_pce_max = np.zeros(len(modes_list))
    err_train_ae_pce = np.zeros(len(modes_list))
    err_train_ae_pce_max = np.zeros(len(modes_list))

    err_ae_base = np.zeros(len(modes_list))
    err_ae_base_max = np.zeros(len(modes_list))

    print("Starting loop over modes")

    for (i,modes) in enumerate(modes_list):
        print(modes)

        pce_model,pce_model_normal,errors,errors_normal,errors_pod_test,active_coeffs,active_coeffs_normal, basis_data, time_POD = gepce(modes,order,partition,parameter_names,regression_type=regression_type,filename="",is_augmented=False,use_saved_basis=False,verbosity=False,angle_mode=angle_mode,num_ds=num_ds)

        checkpoint_path = (f"models/latent_sweep/ae_p5_r{modes}_h256_128_silu_seed42.pt")
        pce_model,pce_model_normal,errors_AE,errors_normal_AE,errors_AE_test_base,active_coeffs,active_coeffs_normal, basis_data_AE, time_AE = gepce_AE(modes,order,partition,parameter_names,regression_type,verbosity=False, angle_mode=angle_mode, checkpoint_path=checkpoint_path, num_ds=num_ds)

        err[i] = errors["test"].mean()
        err_normal[i] = errors_normal["test"].mean()
        err_train[i] = errors["train"].mean()
        err_train_normal[i] = errors_normal["train"].mean()

        err_max[i] = errors["test"].max()
        err_normal_max[i] = errors_normal["test"].max()
        err_train_max[i] = errors["train"].max()
        err_train_normal_max[i] = errors_normal["train"].max()

        err_pod[i] = errors_pod_test.mean()
        err_pod_max[i] = errors_pod_test.max()

        err_ae[i] = errors_AE["test"].mean()
        err_ae_max[i] = errors_AE["test"].max()
        err_train_ae[i] = errors_AE["train"].mean()
        err_train_ae_max[i] = errors_AE["train"].max()

        err_ae_pce[i] = errors_normal_AE["test"].mean()
        err_ae_pce_max[i] = errors_normal_AE["test"].max()
        err_train_ae_pce[i] = errors_normal_AE["train"].mean()
        err_train_ae_pce_max[i] = errors_normal_AE["train"].max()

        err_ae_base[i] = errors_AE_test_base.mean()
        err_ae_base_max[i] = errors_AE_test_base.max()

        times_POD[i] = time_POD
        times_AE[i] = time_AE

    fig, axes = plt.subplots(nrows=1,ncols=3,figsize=(12, 5))

    # Left subplot
    axes[0].plot(modes_list, err, "+-", label="POD-GEPCE")
    axes[0].plot(modes_list, err_normal, "*-", label="POD-PCE")
    axes[0].plot(modes_list, err_pod, "o-", label="POD base")
    axes[0].plot(modes_list, err_ae, "+-", label="AE-GEPCE")
    axes[0].plot(modes_list, err_ae_pce, "*-", label="AE-PCE")
    axes[0].plot(modes_list, err_ae_base, "o-", label="AE base")
    axes[0].axhline(y=err_no_POD, color="red", linestyle="--", label="ref GEPCE")
    axes[0].axhline(y=err_normal_no_POD, color="black", linestyle="--", label="ref PCE")

    axes[0].set_xlabel("Number of modes")
    axes[0].set_ylabel("Relative testing error")
    axes[0].set_title("Mean Approximation error")
    axes[0].grid(True)
    axes[0].legend()
    
    axes[1].plot(modes_list, err_max, "+-", label="POD-GEPCE")
    axes[1].plot(modes_list, err_normal_max, "*-", label="POD-PCE")
    axes[1].plot(modes_list, err_pod_max, "o-", label="POD base")
    axes[1].plot(modes_list, err_ae_max, "+-", label="AE-GEPCE")
    axes[1].plot(modes_list, err_ae_pce_max, "*-", label="AE-PCE")
    axes[1].plot(modes_list, err_ae_base_max, "o-", label="AE base")
    axes[1].axhline(y=err_no_POD_max, color="red", linestyle="--", label="ref GEPCE")
    axes[1].axhline(y=err_normal_no_POD_max, color="black", linestyle="--", label="ref PCE")

    axes[1].set_xlabel("Number of modes")
    axes[1].set_ylabel("Relative testing error")
    axes[1].set_title("Max Approximation error")
    axes[1].grid(True)
    axes[1].legend()

    axes[2].plot(modes_list, times_POD, "+-", label="POD")
    axes[2].plot(modes_list, times_AE, "+-", label="AE")
    axes[2].axhline(y=time_taken_no_POD, color="red", linestyle="--", label="ref")
    axes[2].set_xlabel("Number of modes")
    axes[2].set_ylabel("Time [s]")
    axes[2].set_title("Computation time")
    axes[2].grid(True)
    axes[2].legend()

    fig.suptitle(f"Order{order}, Parameter names: {', '.join(parameter_names)}, training samples {n_train}, testing samples {n_test}")
    fig.tight_layout()

    results_directory = Path("data/ae_pod_comparison_results")
    results_directory.mkdir(parents=True, exist_ok=True)

    plot_data = pd.DataFrame({
        "modes": modes_list,
        # Mean testing errors
        "pod_gepce_mean_error": err,
        "pod_pce_mean_error": err_normal,
        "pod_base_mean_error": err_pod,
        "ae_gepce_mean_error": err_ae,
        "ae_pce_mean_error": err_ae_pce,
        "ae_base_mean_error": err_ae_base,
        "no_compression_gepce_mean_error": err_no_POD,
        "no_compression_pce_mean_error": err_normal_no_POD,
        # Maximum testing errors
        "pod_gepce_max_error": err_max,
        "pod_pce_max_error": err_normal_max,
        "pod_base_max_error": err_pod_max,
        "ae_gepce_max_error": err_ae_max,
        "ae_pce_max_error": err_ae_pce_max,
        "ae_base_max_error": err_ae_base_max,
        "no_compression_gepce_max_error": err_no_POD_max,
        "no_compression_pce_max_error": err_normal_no_POD_max,
        # Computation times
        "pod_time_s": times_POD,
        "ae_time_s": times_AE,
        "no_compression_time_s": time_taken_no_POD,
        # Experiment metadata
        "order": order,
        "partition": partition,
        "num_train": n_train,
        "num_test": n_test,
        "parameters": ", ".join(parameter_names),
    })

    csv_path = results_directory / f"pod_ae_comparison_p{dataset_num_params}_order{order}_part{partition}.csv"
    plot_data.to_csv(csv_path, index=False)
    print(f"Saved plotting data to {csv_path}")

    plt.show()

def compare_timings(parameter_names, dataset_num_params, num_ds, regression_type, modes, orders, partitions, n_runs, angle_mode=SAMPLED_ANGLE_MODE):
    from pce import normalize_design_columns, create_regression_model, calculate_gamma
    from random import randrange

    def compute_training_time_pce(order, basis_data):
        # unpack basis data
        x_train = basis_data["x_train"]
        y_train_red = basis_data["y_train_red"]

        start = time.time()
        pce_basis, Psi, Psi_der = build_design_matrices(order, x_train, with_derivatives=False)
        
        n_coefficients = Psi.shape[1]
        n_modes = y_train_red.shape[1]

        # init coeff matrices
        pce_coeffs = np.zeros((n_coefficients, n_modes))
        Psi_scaled, pce_column_scales = normalize_design_columns(Psi)

        for mode in range(n_modes):
            pce_regressor = create_regression_model(regression_type)
            pce_regressor.fit(Psi_scaled,y_train_red[:, mode])

            # Transform coefficients back to the original polynomial basis.
            pce_coeffs[:, mode] = pce_regressor.coef_ / pce_column_scales

        pce_model = cp.sum(pce_basis[:, None] * pce_coeffs,axis=0)
             
        end = time.time()
        return end - start, pce_model

    def compute_training_time_gepce(order, basis_data):

        # unpack basis data
        x_train = basis_data["x_train"]
        y_train_red = basis_data["y_train_red"]
        y_train_der_red = basis_data["y_train_der_red"]
        
        start = time.time()
        pce_basis, Psi, Psi_der = build_design_matrices(order, x_train)
        
        n_coefficients = Psi.shape[1]
        n_modes = y_train_red.shape[1]

        # init coeff matrices
        gepce_coeffs = np.zeros((n_coefficients, n_modes))
        gamma_values = np.zeros(n_modes)

        for mode in range(n_modes):
            # compute gamme for every mode
            gamma = calculate_gamma(y_train_red[:, mode],y_train_der_red[:, mode])
            gamma_values[mode] = gamma

            Psi_aug = np.vstack([Psi,gamma * Psi_der,])
            y_aug = np.concatenate([y_train_red[:, mode],gamma * y_train_der_red[:, mode]])

            # column normalization
            Psi_aug_scaled, gepce_column_scales = normalize_design_columns(Psi_aug)

            gepce_regressor = create_regression_model(regression_type)
            gepce_regressor.fit(Psi_aug_scaled, y_aug)

            # Transform coefficients back to the original polynomial basis.
            gepce_coeffs[:, mode] = gepce_regressor.coef_ / gepce_column_scales

        gepce_model = cp.sum(pce_basis[:, None] * gepce_coeffs,axis=0)

        end = time.time()
        return end - start, gepce_model

    def compute_prediction_time(model, basis_data):
        # extract data
        x_test = basis_data["x_test"]
        # choose random test sample
        rnd_index = randrange(x_test.shape[0])    
        x_one = x_test[rnd_index:rnd_index + 1, :]

        start = time.time()
        y = predict_pod_fields(model, x_one, basis_data["U"], basis_data["feature_means"], basis_data["feature_stds"])
        end = time.time()
        return end - start

    training_times_pce = np.zeros((len(orders), len(partitions), n_runs))
    prediction_times_pce = np.zeros((len(orders), len(partitions), n_runs))
    training_times_gepce = np.zeros((len(orders), len(partitions), n_runs))
    prediction_times_gepce = np.zeros((len(orders), len(partitions), n_runs))
    samples = np.zeros((len(orders), len(partitions)))

    for i, order in enumerate(orders):
        for j, partition in enumerate(partitions):
            # load data
            filename = (f"data/basis_p{dataset_num_params}_part{partition}_modes{modes}_aug{int(False)}.mat")
            basis_data, time_basis = get_basis(filename, partition, parameter_names, modes, use_saved_basis=True, is_augmented=False, dataset_num_params=dataset_num_params, num_ds=num_ds, angle_mode=angle_mode)
            
            samples[i,j] = basis_data["x_train"].shape[0]

            # warmup computations pce
            warmup_pce, pce_model = compute_training_time_pce(order, basis_data)
            warmup_pred_pce = compute_prediction_time(pce_model, basis_data)

            for n in range(n_runs):
                training_times_pce[i,j,n], _ = compute_training_time_pce(order, basis_data)
                prediction_times_pce[i,j,n] = compute_prediction_time(pce_model, basis_data)

            # warmup computations pce
            warmup_gepce, gepce_model = compute_training_time_gepce(order, basis_data)
            warmup_pred_gepce = compute_prediction_time(gepce_model, basis_data)

            for n in range(n_runs):
                training_times_gepce[i,j,n], _ = compute_training_time_gepce(order, basis_data)
                prediction_times_gepce[i,j,n] = compute_prediction_time(gepce_model, basis_data)

    # compute average values
    trainig_times_pce_averaged = np.average(training_times_pce, axis=2)
    trainig_times_gepce_averaged = np.average(training_times_gepce, axis=2)

    prediction_times_pce_averaged = np.average(prediction_times_pce, axis=2)
    prediction_times_gepce_averaged = np.average(prediction_times_gepce, axis=2)

    return trainig_times_pce_averaged, prediction_times_pce_averaged, trainig_times_gepce_averaged, prediction_times_gepce_averaged, samples

def compare_training_and_prediction_times(n_runs=10):
    import csv
    # run for small dataset
    parameter_names = ["alpha","MH","MW","MAG","Theta1"]
    dataset_num_params = 5
    num_ds = 623
    regression_type = 'lasso'
    modes = 32
    orders = np.array([2,4,6])
    partitions = np.array([0,1,2,3])

    training_pce, prediction_pce, training_gepce, prediction_gepce, samples = compare_timings(parameter_names, dataset_num_params, num_ds, regression_type, modes, orders, partitions, n_runs, angle_mode=SAMPLED_ANGLE_MODE)

    output_path = Path(f"data/timings_p{dataset_num_params}.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow([
            "order", "partition", "samples",
            "pce_training_mean_s", "gepce_training_mean_s",
            "training_ratio", "pce_prediction_mean_s",
            "gepce_prediction_mean_s",
        ])

        for i, order in enumerate(orders):
            for j, partition in enumerate(partitions):
                pce_time = float(training_pce[i, j])
                gepce_time = float(training_gepce[i, j])
                writer.writerow([
                    int(order),
                    int(partition),
                    int(samples[i, j]),
                    pce_time,
                    gepce_time,
                    gepce_time / pce_time,
                    float(prediction_pce[i, j]),
                    float(prediction_gepce[i, j]),
                ])

    # run for big dataset
    parameter_names = ALL_PARAMETERS
    dataset_num_params = 10
    num_ds = 2978
    regression_type = 'lasso'
    modes = 32
    orders = np.array([2,4])
    partitions = np.array([0,1,2,3])
    
    training_pce, prediction_pce, training_gepce, prediction_gepce, samples = compare_timings(parameter_names, dataset_num_params, num_ds, regression_type, modes, orders, partitions, n_runs, angle_mode=SAMPLED_ANGLE_MODE)
    
    output_path = Path(f"data/timings_p{dataset_num_params}.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with output_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow([
                "order", "partition", "samples",
                "pce_training_mean_s", "gepce_training_mean_s",
                "training_ratio", "pce_prediction_mean_s",
                "gepce_prediction_mean_s",
            ])
    
        for i, order in enumerate(orders):
            for j, partition in enumerate(partitions):
                pce_time = float(training_pce[i, j])
                gepce_time = float(training_gepce[i, j])
                writer.writerow([
                    int(order),
                    int(partition),
                    int(samples[i, j]),
                    pce_time,
                    gepce_time,
                    gepce_time / pce_time,
                    float(prediction_pce[i, j]),
                    float(prediction_gepce[i, j]),
                ])

def run_all_experiments():
    # comparison of reconstruction errors POD vs AE, Fig. 4.1
    compare_num_modes_POD_AE()

    # gridded angle case, Fig. 4.2
    compare_gepce(orders=np.array([6]), modes=14, parameter_names=["MW","MH","MAG","Theta1"], partitions=[0,1,2,3], num_ds=163, regression_type="lasso", angle_mode=LEGACY_ANGLE_MODE)

    # POD-(GE)-PCE parameter study on D5, Fig. 4.3, Tab. 4.3
    orders = np.array([2,4,6])
    modes_array = np.array([16,32,48])
    regression_types = ["ols", "lasso"]
    for i, modes in enumerate(modes_array):
        for j, regression_type in enumerate(regression_types):
            compare_gepce(orders, modes, parameter_names=["alpha","MW","MH","MAG","Theta1"], partitions=[0,1,2,3], num_ds=623, regression_type=regression_type, angle_mode=SAMPLED_ANGLE_MODE)

    # comparison of AE-(GE)-PCE and POD-(GE)-PCE and no compression, Fig. 4.4
    compare_POD_AE(order=3)

    # Fig. 4.5: run fieldPlot.m for sample 572

    # ten-parameter case, Fig. 4.6
    compare_gepce(orders=np.array([5]), modes=32, parameter_names=ALL_PARAMETERS, partitions=[0,1,2,3,4], num_ds=2978, regression_type="lasso", angle_mode=SAMPLED_ANGLE_MODE)

    # comprison of training and preediction times, Tab. 4.4, Tab. D.1
    compare_training_and_prediction_times(n_runs=10)

    # Fig. 4.7: run conformalPrediction.py

def main():
    #run_all_experiments()
    compare_training_and_prediction_times(n_runs=10)

if __name__ == "__main__":
    main()
