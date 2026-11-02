import numpy as np
import scipy.io as spio
import os
import time
from dataLoadingHelper import LEGACY_ANGLE_MODE, load_data

POD_BASIS_KEYS = (
    "x_train",
    "x_val",
    "x_test",
    "Y_train",
    "y_train_red",
    "y_train_der_red",
    "Y_val",
    "Y_test",
    "feature_means",
    "feature_stds",
    "U",
    "K_mat",
)

def perform_pod(data_matrix, K, num_modes):
    """ Computes the K-weighted POD basis based on the method of snaphots. Returns the modes as well as the eigenvalues """
    if data_matrix.ndim != 2:
        raise ValueError(f"data_matrix must have shape state_dimension x n_snapshots, but has shape {data_matrix.shape}.")

    # check if number of modes is valid
    if not 1 <= num_modes <= data_matrix.shape[1]:
        raise ValueError(f"num_modes={num_modes} must be between 1 and the number of snapshots ({data_matrix.shape[1]}).")

    # 1. Weighted correlation matrix in snapshot space
    C = data_matrix.T @ K @ data_matrix  # (n_snapshots, n_snapshots)

    # 2. Solve eigenproblem
    eigvals, eigvecs = np.linalg.eigh(C)

    # 3. Sort eigenvalues descending
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # 4. Project back to coefficient space
    modes = data_matrix @ eigvecs  # (n_coeffs, n_snapshots)
    if modes.shape[1] < num_modes:
        raise ValueError(f"Number of training samples {modes.shape[1]} is smaller than number of requested modes {num_modes}")
    
    # 5. Normalize each mode
    for i in range(num_modes):
        phi = modes[:, i]
        norm = np.sqrt(eigvals[i])
        if norm > 1e-12:
            modes[:, i] = phi / norm
        else:
            modes[:, i] = 0.0

    return modes[:,:num_modes], eigvals[:num_modes]

# Project data onto the top k POD modes
def project_onto_pod_modes(data_matrix, U, K):
    """Projects some data_matrix onto the K-weighted POD basis."""
    return U.T @ K @ data_matrix

def project_derivatives_onto_pod(derivatives, U, K_mat):
    """ Project derivatives onto the POD basis.

    derivatives :  (state_dimension, n_parameters, n_samples)
    U : ndarray
    K_mat : ndarray or sparse matrix

    Shape: (n_parameters * n_samples, n_modes).
        The ordering is: parameter 1, all samples; parameter 2, all samples; ...
    """
    # dimension checks
    if derivatives.ndim != 3:
        raise ValueError(f"derivatives must have shape state_dimension x n_parameters x n_samples, but has shape {derivatives.shape}.")
    if derivatives.shape[0] != U.shape[0]:
        raise ValueError(f"Derivative state dimension {derivatives.shape[0]} does not match the POD basis dimension {U.shape[0]}.")

    n_parameters = derivatives.shape[1]
    projected_derivatives = [project_onto_pod_modes(derivatives[:, i_parameter, :], U, K_mat,).T for i_parameter in range(n_parameters)]

    # Initial shape:
    # (n_parameters, n_samples, n_modes)
    projected_derivatives = np.stack(projected_derivatives,axis=0)

    n_modes = projected_derivatives.shape[-1]

    return projected_derivatives.reshape(-1, n_modes)

# load saved basis
def load_saved_basis(filename):
    """Load previously computed POD basis and reduced training data."""
    data = spio.loadmat(filename)
    missing_keys = [key for key in POD_BASIS_KEYS if key not in data]
    if missing_keys:
        raise KeyError(f"Saved basis file is missing entries: {missing_keys}")

    return {key: data[key] for key in POD_BASIS_KEYS}

# save basis data
def save_basis(filename, basis_data):
    """Save POD and reduced training data."""
    if not filename:
        return
    directory = os.path.dirname(filename)
    if directory:
        os.makedirs(directory, exist_ok=True)
    missing_keys = [key for key in POD_BASIS_KEYS if key not in basis_data]
    if missing_keys:
        raise KeyError(f"Basis data is missing entries: {missing_keys}")
    spio.savemat(filename, {key: basis_data[key] for key in POD_BASIS_KEYS})

def get_basis(filename, partition, parameter_names, modes, use_saved_basis, is_augmented, dataset_num_params=None, *, num_ds=163, augmentation_step=1e-3, angle_mode=LEGACY_ANGLE_MODE,dataset_csv=None):

    if use_saved_basis:
        if not filename:
            raise ValueError("A filename is required when use_saved_basis=True.")
        if os.path.exists(filename):
            print(f"Loading saved basis from {filename}")
            return load_saved_basis(filename), 0.0
         
    # load data from files
    x_train, x_val, x_test, data_train, ders_train, data_val, data_test, K_mat, ids_train, ids_val, ids_test  = load_data(
            num_ds=num_ds, partition=partition, parameter_names=parameter_names, dataset_num_params=dataset_num_params, angle_mode=angle_mode,dataset_csv=dataset_csv)

    start = time.time()
    # augment training data with directional derivatives if desired
    if is_augmented:
        pod_snapshots  = data_augmentation(data_train, ders_train, augmentation_step)
    else:
        pod_snapshots  = data_train
            
    # compute POD modes for the original and augmented datasets
    U, eigenvalues = perform_pod(pod_snapshots, K_mat, modes)
    max_imaginary_part = np.max(np.abs(np.imag(U)))
    if max_imaginary_part > 1e-12:
        print(f"Warning: POD modes contain imaginary parts. Maximum magnitude: {max_imaginary_part:.3e}")
    U = np.real_if_close(U, tol=1000)
    if np.iscomplexobj(U):
        raise ValueError("POD modes have significant imaginary components.")
    U = np.asarray(U, dtype=float)

    # project data onto POD modes
    data_train_red = project_onto_pod_modes(data_train, U, K_mat).T
    # Project parameter derivatives
    y_train_der_red = project_derivatives_onto_pod(ders_train,U,K_mat,)

    # Normalize POD coefficients
    feature_means = np.mean(data_train_red,axis=0,keepdims=True)
    feature_stds = np.std(data_train_red,axis=0,keepdims=True)

    y_train_red = (data_train_red - feature_means) / feature_stds
    # Derivatives are not centered because the derivative of the mean is zero
    y_train_der_red = (y_train_der_red / feature_stds)

    basis_data = {
        "x_train": x_train,
        "x_val": x_val,
        "x_test": x_test,
        "Y_train": data_train,
        "y_train_red": y_train_red,
        "y_train_der_red": y_train_der_red,
        "Y_val": data_val,
        "Y_test": data_test,
        "feature_means": feature_means,
        "feature_stds": feature_stds,
        "U": U,
        "K_mat": K_mat,
    }

    end = time.time()
    time_taken = end - start
    save_basis(filename, basis_data)

    return basis_data, time_taken

def data_augmentation(data_train, ders_train, h=1e-3):
    data_train_aug = data_train.copy()
    snaps_der_train = ders_train.reshape(ders_train.shape[0], -1)

    # augment training data with directional derivatives
    for ii in range(data_train.shape[1]):
        base = data_train[:, ii:ii+1]                # (ndof, 1)
        grads = ders_train[:, :, ii]                 # (ndof, n_params)
        vect_add_pos = base + h * grads
        vect_add_neg = base - h * grads
        data_train_aug = np.concatenate((data_train_aug, vect_add_pos, vect_add_neg), axis=1)

    return data_train_aug
