import numpy as np
import os
import pandas as pd
import scipy.io as spio
import mat73

from parameters import PARAMETER_RANGES

LEGACY_ANGLE_MODE = "legacy_grid"
SAMPLED_ANGLE_MODE = "sampled_alpha"

def matlab_strings(values):
    """Convert common MATLAB string/cell layouts to a flat Python list."""
    if isinstance(values, str):
        return [values]

    names = []
    for value in np.asarray(values, dtype=object).ravel():
        array = np.asarray(value).squeeze()

        if array.ndim == 0:
            item = array.item()
            if isinstance(item, bytes):
                item = item.decode()
            names.append(str(item))
        elif array.dtype.kind in {"U", "S"}:
            names.append("".join(array.astype(str).ravel()))
        else:
            names.append(str(array))

    return names

def normalize_data(data, ranges):
    """Normalizes data given their ranges."""
    normalized = np.zeros_like(data)
    for i, (min_val, max_val) in enumerate(ranges):
        normalized[i, :] = (data[i, :] - min_val) / (max_val - min_val)
    return normalized

def load_mat_file(filename):
    """Load a complete MATLAB file, including v7.3 files."""
    try:
        # MATLAB files older than v7.3
        return spio.loadmat(filename)
    except NotImplementedError:
        # MATLAB v7.3 files
        return mat73.loadmat(filename)

def _load_sample_table(num_ds, partition, parameter_names, dataset_num_params, *, dataset_csv=None):
    """Read the CSV, select one Sobol partition, and normalize its inputs."""
    # checks for unknown parameters or wrong names
    parameter_names = list(parameter_names)
    if "RotationAngle" in parameter_names:
        raise ValueError("Use 'alpha', not the internal name 'RotationAngle'.")
    if not parameter_names or len(parameter_names) != len(set(parameter_names)):
        raise ValueError("parameter_names must contain unique parameter names.")
    unknown = set(parameter_names) - set(PARAMETER_RANGES)
    if unknown:
        raise ValueError(f"Unknown parameters: {sorted(unknown)}")
    if dataset_num_params is None:
        dataset_num_params = len(parameter_names)

    # build filenames
    stem = f"p_{dataset_num_params}_samp_{num_ds}"
    # this is annyoing but I need this to allow for my CP evaluation data with seeds in the name
    if dataset_csv is None:
        csv_file = os.path.join("data", f"dataset_{stem}.csv")
    else:
        csv_file = os.fspath(dataset_csv)

    csv_dir = os.path.dirname(csv_file)
    csv_stem = os.path.splitext(os.path.basename(csv_file))[0]

    mat_file = os.path.join(csv_dir, f"{csv_stem}.mat")
    matrix_stem = csv_stem.replace("dataset_", "error_matrices_", 1)
    matrix_file = os.path.join(csv_dir, f"{matrix_stem}.mat")

    table = pd.read_csv(csv_file)
    required_columns = {"dataset", "sobol_id", "global_id", *parameter_names}
    missing = required_columns - set(table.columns)
    if missing:
        raise ValueError(f"Missing CSV columns: {sorted(missing)}")

    # Check whether all parameters are contained in the allowed ranges. This also catches accidentally storing alpha in degrees instead of radians.
    ranges = [PARAMETER_RANGES[name] for name in parameter_names]
    for name, (lower, upper) in zip(parameter_names, ranges):
        if not table[name].between(lower, upper).all():
            raise ValueError(f"Values of '{name}' are outside [{lower}, {upper}].")

    # check whether the requested sobol partion exists
    table = table.loc[table["sobol_id"] <= partition].copy()
    if table.empty:
        raise ValueError(f"Sobol partition {partition} contains no samples.")

    # MATLAB global_id starts at one, while Python indexing starts at zero.
    global_ids = table["global_id"].to_numpy(dtype=int) - 1
    if np.any(global_ids < 0):
        raise ValueError("global_id must start at one in the CSV file.")

    labels = table["dataset"].to_numpy()
    local_ids = np.arange(len(table))
    train_ids = local_ids[labels == 0]
    val_ids = local_ids[labels == 1]
    test_ids = local_ids[labels == 2]

    parameter_values = table[parameter_names].to_numpy(dtype=float).T
    x = normalize_data(parameter_values, ranges).T
    widths = np.array([upper - lower for lower, upper in ranges])

    return {
        "parameter_names": parameter_names,
        "table": table,
        "mat_file": mat_file,
        "matrix_file": matrix_file,
        "global_ids": global_ids,
        "global_train_ids": global_ids[train_ids],
        "train_ids": train_ids,
        "val_ids": val_ids,
        "test_ids": test_ids,
        "x_train": x[train_ids],
        "x_val": x[val_ids],
        "x_test": x[test_ids],
        "widths": widths,
    }

def load_legacy_angle_grid_data(num_ds, partition, parameter_names, dataset_num_params=None,*,dataset_csv=None):
    """Load the old format: one CSV sample contains all rotor angles."""

    info = _load_sample_table(num_ds, partition, parameter_names, dataset_num_params,dataset_csv=dataset_csv)
    if "alpha" in info["table"].columns:
        raise ValueError("The legacy CSV must not contain an 'alpha' column.")

    mat_data = load_mat_file(info["mat_file"])
    state_store = np.asarray(mat_data["state_store"])
    state_der_store = np.asarray(mat_data["state_der_store"])
    derivative_names = matlab_strings(mat_data["derivative_names"])

    # Fixed legacy layout: angles x dofs x samples [x parameters].
    if state_store.ndim != 3 or state_der_store.ndim != 4:
        raise ValueError("Legacy arrays must have shapes (angles, dofs, samples) and (angles, dofs, samples, parameters).")
    if state_der_store.shape[:3] != state_store.shape:
        raise ValueError("state_store and state_der_store dimensions do not match.")
    if state_der_store.shape[3] != len(derivative_names):
        raise ValueError("derivative_names does not match state_der_store.")
    if info["global_ids"].max() >= state_store.shape[2]:
        raise ValueError("A global_id exceeds the number of stored snapshots.")

    missing = set(info["parameter_names"]) - set(derivative_names)
    if missing:
        raise ValueError(f"Missing derivatives: {sorted(missing)}")

    snapshots = np.take(state_store, info["global_ids"], axis=2)

    # Select one derivative at a time to avoid a large temporary array.
    n_angles, n_dofs = state_store.shape[:2]
    derivatives = np.empty((n_angles, n_dofs, len(info["parameter_names"]), len(info["train_ids"])), dtype=state_der_store.dtype)
    for j, name in enumerate(info["parameter_names"]):
        index = derivative_names.index(name)
        derivatives[:, :, j, :] = np.take(state_der_store[:, :, :, index], info["global_train_ids"], axis=2)

    # All angle fields form one large output vector in the legacy workflow.
    snapshots = snapshots.reshape(n_angles * n_dofs, len(info["global_ids"]))
    derivatives = derivatives.reshape(n_angles * n_dofs, len(info["parameter_names"]), len(info["train_ids"]))

    K_mat = load_mat_file(info["matrix_file"])["Kred_big"]

    # derivative normalization: multiply by the width of each parameter range
    derivatives *= info["widths"][None, :, None]

    return (
        info["x_train"],
        info["x_val"],
        info["x_test"],
        snapshots[:, info["train_ids"]],
        derivatives,
        snapshots[:, info["val_ids"]],
        snapshots[:, info["test_ids"]],
        K_mat,
        info["train_ids"],
        info["val_ids"],
        info["test_ids"],
    )

def load_sampled_alpha_data(num_ds, partition, parameter_names, dataset_num_params=None,*,dataset_csv=None):
    """Load the new format: every CSV row contains one ``alpha`` value."""

    info = _load_sample_table(num_ds, partition, parameter_names, dataset_num_params,dataset_csv=dataset_csv)
    # check if alpha is contained
    if "alpha" not in info["table"].columns:
        raise ValueError("The sampled-angle CSV must contain 'alpha'.")
    if "alpha" not in info["parameter_names"]:
        raise ValueError("The sampled-angle PCE must use 'alpha'.")

    # load stuff
    mat_data = load_mat_file(info["mat_file"])
    state_store = np.asarray(mat_data["state_store"]).squeeze()
    state_der_store = np.asarray(mat_data["state_der_store"]).squeeze()
    derivative_names = matlab_strings(mat_data["derivative_names"])

    # Fixed sampled-alpha layout after removing MATLAB singleton dimensions:
    # dofs x samples [x parameters].
    if state_der_store.ndim == 2 and len(derivative_names) == 1:
        state_der_store = state_der_store[:, :, None]
    if state_store.ndim != 2 or state_der_store.ndim != 3:
        raise ValueError("Sampled-alpha arrays must have shapes (dofs, samples) and (dofs, samples, parameters) after squeezing.")
    if state_der_store.shape[:2] != state_store.shape:
        raise ValueError("state_store and state_der_store dimensions do not match.")
    if state_der_store.shape[2] != len(derivative_names):
        raise ValueError("derivative_names does not match state_der_store.")
    if info["global_ids"].max() >= state_store.shape[1]:
        raise ValueError("A global_id exceeds the number of stored snapshots.")

    missing = set(info["parameter_names"]) - set(derivative_names)
    if missing:
        raise ValueError(f"Missing derivatives: {sorted(missing)}. The angle derivative must be stored under the exact name 'alpha'.")

    snapshots = np.take(state_store, info["global_ids"], axis=1)

    derivatives = np.empty((state_store.shape[0], len(info["parameter_names"]), len(info["train_ids"])), dtype=state_der_store.dtype)
    for j, name in enumerate(info["parameter_names"]):
        index = derivative_names.index(name)
        derivatives[:, j, :] = np.take(state_der_store[:, :, index], info["global_train_ids"], axis=1)

    K_mat = load_mat_file(info["matrix_file"])["Kred"]

    # derivative normalization: multiply by the width of each parameter range
    derivatives *= info["widths"][None, :, None]

    return (
        info["x_train"],
        info["x_val"],
        info["x_test"],
        snapshots[:, info["train_ids"]],
        derivatives,
        snapshots[:, info["val_ids"]],
        snapshots[:, info["test_ids"]],
        K_mat,
        info["train_ids"],
        info["val_ids"],
        info["test_ids"],
    )

def load_data(num_ds,  partition,  parameter_names, dataset_num_params=None, *, angle_mode=LEGACY_ANGLE_MODE,dataset_csv=None):
    """Dispatch explicitly to one of the two supported data representations."""

    if angle_mode == LEGACY_ANGLE_MODE:
        loader = load_legacy_angle_grid_data
    elif angle_mode == SAMPLED_ANGLE_MODE:
        loader = load_sampled_alpha_data
    else:
        raise ValueError(f"Unknown angle_mode '{angle_mode}'. Choose '{LEGACY_ANGLE_MODE}' or '{SAMPLED_ANGLE_MODE}'.")

    return loader(num_ds, partition, parameter_names, dataset_num_params=dataset_num_params,dataset_csv=dataset_csv)
