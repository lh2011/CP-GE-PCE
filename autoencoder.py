import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from pathlib import Path
import optuna
import time

from dataLoadingHelper import LEGACY_ANGLE_MODE, load_data, load_sampled_alpha_data

ARCHITECTURE_OPTIONS = {
    "small": (128, 64),
    "medium": (256, 64),
    "wide": (256, 128),
    #"very-wide": (512, 128),
    #"ultra-wide": (512, 256),
}

class FieldAutoencoder(nn.Module):
    def __init__(self,num_dofs,latent_dim,activation="gelu",hidden_dim_1=128,hidden_dim_2=64):

        super().__init__()

        if latent_dim <= 0:
            raise ValueError("latent_dim must be positive.")

        if latent_dim > hidden_dim_2:
            raise ValueError("latent_dim should not exceed hidden_dim_2, because it would no longer be the narrowest layer.")

        self.num_dofs = num_dofs
        self.latent_dim = latent_dim
        if activation == "gelu":
            activation_fn = nn.GELU()
        elif activation == "relu":
            activation_fn = nn.ReLU()
        elif activation == "silu":
            activation_fn = nn.SiLU()
        elif activation == "tanh":
            activation_fn = nn.Tanh()
        elif activation == "sigm":
            activation_fn = nn.Sigmoid()
        else:
            raise ValueError(f"Unknown activation: {activation}")

        self.encoder = nn.Sequential(
            nn.Linear(num_dofs, hidden_dim_1),
            activation_fn,
            nn.Linear(hidden_dim_1, hidden_dim_2),
            activation_fn,
            nn.Linear(hidden_dim_2, latent_dim),
        )

        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim_2),
            activation_fn,            
            nn.Linear(hidden_dim_2, hidden_dim_1),
            activation_fn,
            nn.Linear(hidden_dim_1, num_dofs),
        )

    def encode(self, fields):
        """Map normalized fields to latent coordinates."""
        return self.encoder(fields)

    def decode(self, latent_variables):
        """Map latent coordinates to normalized fields."""
        return self.decoder(latent_variables)

    def forward(self, fields):
        latent_variables = self.encode(fields)
        reconstructed_fields = self.decode(latent_variables)

        return reconstructed_fields

def scipy_to_torch_sparse(matrix, device):
    """Converts a scipy matrix too a COO torch matrix."""
    matrix = matrix.tocoo()
    indices = torch.tensor(np.vstack((matrix.row, matrix.col)), dtype=torch.long)
    values = torch.tensor(matrix.data, dtype=torch.float32)
    torch_matrix = torch.sparse_coo_tensor(indices=indices, values=values, size=matrix.shape, dtype=torch.float32)
    return torch_matrix.coalesce().to(device)

def save_autoencoder_checkpoint(model, checkpoint_path, mean_field, field_scale, selected_ids, val_ids, best_epoch, validation_error, seed, hyperparameters):
    """Saves a checkpoint (ergo the model) of the autoencoder with all required data."""
    checkpoint_path = Path(checkpoint_path)
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint = {
        "model_state_dict": {
            name: value.detach().cpu() for name, value in model.state_dict().items()
        },
        "model_config": {
            "num_dofs": model.num_dofs,
            "latent_dim": model.latent_dim,
            "hidden_dim_1": hyperparameters["hidden_dim_1"],
            "hidden_dim_2": hyperparameters["hidden_dim_2"],
            "activation": hyperparameters["activation"],
        },
        "preprocessing": {
            "mean_field": torch.as_tensor(mean_field, dtype=torch.float64),
            "field_scale": float(field_scale),
        },
        "training": {
            "seed": int(seed),
            "best_epoch": int(best_epoch),
            "validation_error": float(validation_error),
            "selected_ids": torch.as_tensor(selected_ids, dtype=torch.long),
            "val_ids": torch.as_tensor(val_ids, dtype=torch.long),
            "learning_rate": float(hyperparameters["learning_rate"]),
            "weight_decay": float(hyperparameters["weight_decay"]),
            "batch_size": int(hyperparameters["batch_size"]),
        },
    }
    torch.save(checkpoint, checkpoint_path)
    print(f"Saved model to {checkpoint_path}")

def load_autoencoder_checkpoint(checkpoint_path, device=None):
    """Load a saved checkpoint of the field autoencoder."""
    checkpoint_path = Path(checkpoint_path)
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    model = FieldAutoencoder(**checkpoint["model_config"]).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    mean_field = checkpoint["preprocessing"]["mean_field"].cpu().numpy()
    field_scale = checkpoint["preprocessing"]["field_scale"]

    return {
        "model": model,
        "mean_field": mean_field,
        "field_scale": field_scale,
        "training": checkpoint["training"],
        "model_config": checkpoint["model_config"],
        "device": device,
    }

def get_basis_AE(partition, parameter_names, modes, dataset_num_params=None, *, checkpoint_path, num_ds=163, angle_mode=LEGACY_ANGLE_MODE):
    """Loads the data and then constructs the latent basis representation using the Autoencoder specified in the checkpoint_path."""
    # load data from files
    x_train, x_val, x_test, data_train, ders_train, data_val, data_test, K_mat, ids_train, ids_val, ids_test  = load_data(
            num_ds=num_ds, partition=partition, parameter_names=parameter_names, dataset_num_params=dataset_num_params, angle_mode=angle_mode)

    start = time.time()

    # load Autoencoder
    loaded_ae = load_autoencoder_checkpoint(checkpoint_path)
    model = loaded_ae["model"]
    mean_field = loaded_ae["mean_field"]
    field_scale = loaded_ae["field_scale"]
    device = loaded_ae["device"]

    if model.latent_dim != modes:
        raise ValueError(f"Checkpoint has latent_dim={model.latent_dim}, but modes={modes} was requested.")
    if model.num_dofs != data_train.shape[0]:
        raise ValueError(f"Checkpoint expects {model.num_dofs} DOFs, but the data contains {data_train.shape[0]}.")

    # Normalize fields exactly as during AE training.
    network_train = np.ascontiguousarray(((data_train - mean_field) / field_scale).T, dtype=np.float32)
    network_train_torch = torch.from_numpy(network_train).to(device)

    model.eval()

    # Encode field snapshots: (samples, DOFs) -> (samples, modes)
    with torch.no_grad():
        latent_train = model.encode(network_train_torch).cpu().numpy()

    if latent_train.shape != (x_train.shape[0], modes):
        raise ValueError(f"Unexpected latent shape {latent_train.shape}; expected {(x_train.shape[0], modes)}.")

    # Mean and standard deviation of each latent coordinate.
    latent_means = np.mean(latent_train, axis=0, keepdims=True)
    latent_stds = np.std(latent_train, axis=0, keepdims=True)

    if not np.all(np.isfinite(latent_stds)) or np.any(latent_stds == 0.0):
        raise ValueError("At least one AE latent coordinate has zero or invalid variance.")

    # PCE targets: standardized latent coordinates.
    y_train_red = (latent_train - latent_means) / latent_stds

    # (DOFs, parameters, samples)
    # -> (parameters, samples, DOFs)
    ders_train = ders_train.transpose(1, 2, 0)

    # encode derivative data
    latent_derivatives = []
    for parameter_derivative in ders_train:
        parameter_derivative = np.asarray(parameter_derivative)

        # d[(u - mean_field) / field_scale] / dmu
        normalized_derivative = np.ascontiguousarray(parameter_derivative / field_scale, dtype=np.float32)
        derivative_torch = torch.from_numpy(normalized_derivative).to(device)

        # dz/dmu = J_encoder(u_normalized) @ du_normalized/dmu
        with torch.no_grad():
            _, latent_derivative = torch.func.jvp(model.encode,(network_train_torch,),(derivative_torch,))
        latent_derivatives.append(latent_derivative.detach().cpu().numpy())

    # Shape: (parameters, samples, modes)
    latent_derivatives = np.stack(latent_derivatives, axis=0)

    # Standardizing z also requires standardizing dz/dmu.
    latent_derivatives_normalized = latent_derivatives / latent_stds

    # Parameter-major ordering:
    # parameter 1/all samples, parameter 2/all samples, ...
    y_train_der_red = latent_derivatives_normalized.reshape(len(parameter_names) * x_train.shape[0],modes)

    basis_data = {
    "x_train": x_train,
    "x_val": x_val,
    "x_test": x_test,
    "Y_train": data_train,
    "Y_val": data_val,
    "Y_test": data_test,
    "y_train_red": y_train_red,
    "y_train_der_red": y_train_der_red,
    "feature_means": latent_means,
    "feature_stds": latent_stds,
    "mean_field": mean_field,
    "field_scale": field_scale,
    "model": model,
    "K_mat": K_mat,
    }

    end = time.time()
    time_taken = end - start

    return basis_data, time_taken

def prepare_autoencoder_data(data_train, data_test, num_train=256, num_val=64, seed=42):
    """This splits the data so that there is enough test testing data and also normalizes it."""
    if num_train > data_train.shape[1]:
        raise ValueError(f"Requested {num_train} training snapshots, but only {data_train.shape[1]} are available.")

    rng = np.random.default_rng(seed)
    selected_ids = rng.choice(data_train.shape[1], size=num_train, replace=False)
    snapshots_train = data_train[:, selected_ids]

    # Compute preprocessing quantities from training data only.
    mean_field = snapshots_train.mean(axis=1, keepdims=True)
    centered_train = snapshots_train - mean_field

    # One global scaling factor.
    field_scale = np.sqrt(np.mean(centered_train**2))
    if not np.isfinite(field_scale) or field_scale <= 0.0:
        raise ValueError(f"Invalid field scale: {field_scale}")

    # PyTorch uses shape (number of samples, number of DOFs).
    network_train = np.ascontiguousarray(((snapshots_train - mean_field) / field_scale).T, dtype=np.float32)
    network_test = np.ascontiguousarray(((data_test - mean_field) / field_scale).T, dtype=np.float32)

    ae_data =  {
        "snapshots_train": snapshots_train,
        "snapshots_test": data_test,
        "network_train": network_train,
        "network_test": network_test,
        "mean_field": mean_field,
        "field_scale": field_scale,
        "selected_ids": selected_ids,
    }

    all_ids = np.arange(data_train.shape[1])
    unused_ids = np.setdiff1d(all_ids, ae_data["selected_ids"])
    if num_val > unused_ids.size:
        raise ValueError("Not enough unused samples.")
    
    np.random.default_rng(seed).shuffle(unused_ids)
    val_ids = unused_ids[:num_val]
    snapshots_val = data_train[:, val_ids]
    network_val = np.ascontiguousarray(((snapshots_val - ae_data["mean_field"]) / ae_data["field_scale"] ).T, dtype=np.float32)
    
    ae_data["snapshots_val"] = snapshots_val
    ae_data["network_val"] = network_val
    ae_data["val_ids"] = val_ids
    
    print(f"Network training data: {network_train.shape}")
    print(f"Network testing data:  {network_test.shape}")
    print(f"Network validation data: {network_val.shape}")
    print(f"Field scale:           {field_scale:.6e}")
    print(f"Normalized training RMS: {np.sqrt(np.mean(network_train**2)):.6f}")

    return ae_data

def prepare_training_setup(*, num_ds=623, partition=4, parameter_names=("alpha", "MH", "MW", "MAG", "Theta1"), dataset_num_params=5, num_train=400, num_val=96, seed=42):
    _, _, _, data_train, _, _, data_test, K_mat, _, _, _ = load_sampled_alpha_data(num_ds=num_ds, partition=partition, parameter_names=parameter_names, dataset_num_params=dataset_num_params)
    
    ae_data = prepare_autoencoder_data(data_train=data_train, data_test=data_test, num_train=num_train, num_val=num_val, seed=seed, )

    torch.manual_seed(seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    K_symmetric = 0.5 * (K_mat + K_mat.T)
    K_torch = scipy_to_torch_sparse(K_symmetric, device)

    train_tensor = torch.from_numpy(ae_data["network_train"])
    training_fields = train_tensor.to(device)
    validation_physical = torch.from_numpy(np.ascontiguousarray(ae_data["snapshots_val"].T, dtype=np.float32)).to(device)

    with torch.no_grad():
        K_training = torch.sparse.mm(K_torch, training_fields.T).T
        E_train = (training_fields * K_training).sum(dim=1).mean().item()
        
        K_validation = torch.sparse.mm(K_torch, validation_physical.T).T
        validation_physical_energies = (validation_physical * K_validation).sum(dim=1)

    return {
        "ae_data": ae_data,
        "device": device,
        "K_mat": K_mat,
        "K_torch": K_torch,
        "E_train": E_train,
        "validation_physical_energies": validation_physical_energies,
        "train_dataset": TensorDataset(train_tensor),
        "seed": seed,
        "dataset_num_params": dataset_num_params,
    }

def reconstruct_fields(model, normalized_fields, mean_field, field_scale, device):
    """Reconstruct the fields given the normalized fields"""
    model.eval()
    fields_tensor = torch.from_numpy(normalized_fields).to(device)
    with torch.no_grad():
        reconstruction_normalized = (model(fields_tensor).cpu().numpy())

    # Convert from (samples, DOFs) back to (DOFs, samples).
    reconstruction = mean_field + field_scale * reconstruction_normalized.T

    return reconstruction

def relative_energy_errors(reference_fields, reconstructed_fields, K_mat):
    """Computes the K-weighted loss normlized the reference field energy."""
    residual = reconstructed_fields - reference_fields
    # note for me: this summation over the dofs to give one value for every sample
    numerator = np.sum(residual * (K_mat @ residual), axis=0)
    denominator = np.sum(reference_fields * (K_mat @ reference_fields), axis=0)
    numerator = np.maximum(numerator, 0.0)
    return np.sqrt(numerator / denominator)

def k_weighted_loss(reconstructed_fields, reference_fields, K_torch, E_reference):
    """Computes the K-weighted loss normlized by the reference Energy."""
    residual = reconstructed_fields - reference_fields
    K_residual = torch.sparse.mm(K_torch, residual.T).T
    residual_energies = torch.sum(residual * K_residual,  dim=1,)
    return residual_energies.mean() / E_reference

def mean_relative_energy_error_torch(reconstructed_normalized, reference_normalized, reference_physical_energies, field_scale, K_torch):
    """Computes the mean relative energy error normalied by the reference energies which have been precomputed in the physical frame."""
    # Convert the normalized residual back to physical scaling.
    residual_physical = field_scale * (reconstructed_normalized - reference_normalized)
    K_residual = torch.sparse.mm(K_torch, residual_physical.T).T

    numerator = torch.sum(residual_physical * K_residual, dim=1)
    numerator = torch.clamp(numerator, min=0.0)

    relative_errors = torch.sqrt(numerator / reference_physical_energies)

    return relative_errors.mean()

def train_fieldautoencoder(latent_dim, hidden_dim_1, hidden_dim_2, activation,batch_size, learning_rate, weight_decay, model_subfolder="final"):
    setup = prepare_training_setup()
    ae_data = setup["ae_data"]
    device = setup["device"]
    seed = setup["seed"]
    dataset_num_params = setup["dataset_num_params"]
    K_mat = setup["K_mat"]
    num_dofs = ae_data["network_train"].shape[1]

    # construct Autoencoder + Loader
    model = FieldAutoencoder(
        num_dofs=num_dofs,
        latent_dim=latent_dim,
        activation=activation,
        hidden_dim_1=hidden_dim_1,
        hidden_dim_2=hidden_dim_2
    ).to(device)
    train_loader = DataLoader(setup["train_dataset"], batch_size=batch_size, shuffle=True, generator=torch.Generator().manual_seed(seed))

    # train autoencoder
    best_validation_error, best_epoch = train_autoencoder(
        model=model,
        train_loader=train_loader,
        validation_data=ae_data["network_val"],
        validation_physical_energies=setup["validation_physical_energies"],
        K_torch=setup["K_torch"],
        E_train=setup["E_train"],
        field_scale=ae_data["field_scale"],
        device=device,
        max_epochs=3000,
        learning_rate=learning_rate,
        weight_decay=weight_decay,
        patience=50,
        min_delta=1e-6,
    )

    checkpoint_hyperparameters = {
        "latent_dim": latent_dim,
        "hidden_dim_1": hidden_dim_1,
        "hidden_dim_2": hidden_dim_2,
        "activation": activation,
        "batch_size": batch_size,
        "learning_rate": learning_rate,
        "weight_decay": weight_decay,
    }
    checkpoint_directory = Path("models") / model_subfolder
    checkpoint_path = checkpoint_directory / f"ae_p{dataset_num_params}_r{latent_dim}_h{hidden_dim_1}_{hidden_dim_2}_{activation}_seed{seed}.pt"
    save_autoencoder_checkpoint(
        model=model,
        checkpoint_path=checkpoint_path,
        mean_field=ae_data["mean_field"],
        field_scale=ae_data["field_scale"],
        selected_ids=ae_data["selected_ids"],
        val_ids=ae_data["val_ids"],
        best_epoch=best_epoch,
        validation_error=best_validation_error,
        seed=seed,
        hyperparameters=(checkpoint_hyperparameters),
    )
    
    # reconstruct training and testing fields and compute relative energy errors
    reconstruction_train = reconstruct_fields(model=model, normalized_fields=ae_data["network_train"], mean_field=ae_data["mean_field"], field_scale=ae_data["field_scale"], device=device)
    reconstruction_test = reconstruct_fields(model=model,normalized_fields=ae_data["network_test"],mean_field=ae_data["mean_field"], field_scale=ae_data["field_scale"], device=device)

    errors_train = relative_energy_errors(reference_fields=ae_data["snapshots_train"], reconstructed_fields=reconstruction_train, K_mat=K_mat)
    errors_test = relative_energy_errors(reference_fields=ae_data["snapshots_test"], reconstructed_fields=reconstruction_test, K_mat=K_mat)

    print(f"Mean training error: {100 * errors_train.mean():.4f} %")
    print(f"Mean testing error:  {100 * errors_test.mean():.4f} %")
    print(f"Median testing error: {100 * np.median(errors_test):.4f} %")
    print(f"Maximum testing error: {100 * errors_test.max():.4f} %")

    return model

def train_autoencoder(model, train_loader, validation_data, validation_physical_energies, K_torch, E_train, field_scale, device, trial=None, max_epochs=3000, learning_rate=1e-3, weight_decay=0.0, patience=50, min_delta=1e-6):

    optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    validation_fields = torch.from_numpy(validation_data).to(device)

    best_validation_error = np.inf
    best_epoch = 0
    best_model_state = None
    epochs_without_improvement = 0

    for epoch in range(max_epochs):
        model.train()
        total_training_loss = 0.0
        number_samples = 0

        for (fields,) in train_loader:
            fields = fields.to(device)
            optimizer.zero_grad(set_to_none=True)

            reconstructed_fields = model(fields)
            loss = k_weighted_loss(reconstructed_fields,fields, K_torch, E_train)
            loss.backward()
            optimizer.step()

            batch_size = fields.shape[0]
            total_training_loss += loss.item() * batch_size
            number_samples += batch_size

        mean_training_loss = total_training_loss / number_samples
        
        model.eval()
        with torch.no_grad():
            reconstructed_validation = model(validation_fields)

            validation_error = mean_relative_energy_error_torch(reconstructed_normalized=reconstructed_validation, reference_normalized=validation_fields, reference_physical_energies=validation_physical_energies, field_scale=field_scale,  K_torch=K_torch).item()
            
        # Make intermediate results available to Optuna.
        if trial is not None:
            trial.report(validation_error, epoch)
            if trial.should_prune():
                raise optuna.TrialPruned()

        if (validation_error < best_validation_error - min_delta):
            best_validation_error = validation_error
            best_epoch = epoch + 1
            best_model_state = {name: value.detach().cpu().clone() for name, value in model.state_dict().items()}
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1

        if epoch == 0 or (epoch + 1) % 100 == 0:
            print(f"Epoch {epoch + 1:4d}: train={mean_training_loss:.6e}, validation={100*validation_error:.4f} %")
        if epochs_without_improvement >= patience:
            break

    if best_model_state is None:
        raise RuntimeError( "No valid model state was obtained.")

    model.load_state_dict(best_model_state)

    return best_validation_error, best_epoch

def hyperparameter_optimization():
    setup = prepare_training_setup()
    ae_data = setup["ae_data"]
    device = setup["device"]
    seed = setup["seed"]
    train_dataset = setup["train_dataset"]
    num_dofs = ae_data["network_train"].shape[1]

    def objective(trial):
        # Use the same initialization and shuffling conditions to make trials more directly comparable.
        torch.manual_seed(seed)

        architecture_name = trial.suggest_categorical("architecture",  list(ARCHITECTURE_OPTIONS))
        hidden_dim_1, hidden_dim_2 = (ARCHITECTURE_OPTIONS[architecture_name])
        activation = trial.suggest_categorical("activation", ["gelu", "silu"])
        learning_rate = trial.suggest_float("learning_rate", 1e-5, 3e-3, log=True)
        weight_decay = trial.suggest_float( "weight_decay", 1e-8, 1e-3, log=True)
        batch_size = trial.suggest_categorical("batch_size", [16, 32, 64])

        train_generator = torch.Generator()
        train_generator.manual_seed(seed)
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, generator=train_generator)
        model = FieldAutoencoder(num_dofs=num_dofs, latent_dim=32, activation=activation, hidden_dim_1=hidden_dim_1, hidden_dim_2=hidden_dim_2).to(device)

        best_validation_error, best_epoch = (
            train_autoencoder(
                model=model,
                train_loader=train_loader,
                validation_data=ae_data["network_val"],
                validation_physical_energies=setup["validation_physical_energies"],
                K_torch=setup["K_torch"],
                E_train=setup["E_train"],
                field_scale=ae_data["field_scale"],
                device=device,
                trial=trial,
                max_epochs=3000,
                learning_rate=learning_rate,
                weight_decay=weight_decay,
                patience=50,
                min_delta=1e-6,
            )
        )

        trial.set_user_attr("best_epoch", best_epoch)
        return best_validation_error

    Path("models").mkdir(parents=True, exist_ok=True)

    sampler = optuna.samplers.TPESampler(seed=42, n_startup_trials=10)
    pruner = optuna.pruners.MedianPruner(n_startup_trials=10, n_warmup_steps=50, interval_steps=10)

    study = optuna.create_study(
        study_name="ae_p5_r32",
        storage="sqlite:///models/ae_optuna_p5_r32.db",
        load_if_exists=True,
        direction="minimize",
        sampler=sampler,
        pruner=pruner,
    )

    """
    study = optuna.create_study(
        study_name="ae_p5_r32_v2",
        storage="sqlite:///models/ae_optuna_p5_r32_v2.db",
        load_if_exists=True,
        direction="minimize",
        sampler=sampler,
        pruner=pruner,
    )
    """

    study.optimize(
        objective,
        n_trials=50,
        gc_after_trial=True,
        show_progress_bar=True,
    )

    print("\nBest trial")
    print(f"Validation error: {100*study.best_value:.4f} %")
    print(f"Parameters:       {study.best_params}")
    print(f"Best epoch: {study.best_trial.user_attrs['best_epoch']}")

    best_params = study.best_params
    hidden_dim_1, hidden_dim_2 = ARCHITECTURE_OPTIONS[best_params["architecture"]]
    final_hyperparameters = {
        "hidden_dim_1": hidden_dim_1,
        "hidden_dim_2": hidden_dim_2,
        "activation": best_params["activation"],
        "learning_rate": best_params["learning_rate"],
        "weight_decay": best_params["weight_decay"],
        "batch_size": best_params["batch_size"],
    }

    print(final_hyperparameters)

def train_autoencoders_for_comparison(modes_list):

    # load data from hyperparameter optimization
    database_path = "models/ae_optuna_p5_r32.db"
    study_name = "ae_p5_r32"

    #database_path = "models/ae_optuna_p5_r32_v2.db"
    #study_name = "ae_p5_r32_v2"

    database_path = Path(database_path).resolve()
    storage_url = f"sqlite:///{database_path.as_posix()}"

    study = optuna.load_study(study_name=study_name, storage=storage_url)
    best_params = study.best_params
    architecture_name = best_params["architecture"]
    hidden_dim_1, hidden_dim_2 = (ARCHITECTURE_OPTIONS[architecture_name])

    configuration = {
        "hidden_dim_1": hidden_dim_1,
        "hidden_dim_2": hidden_dim_2,
        "activation": best_params["activation"],
        "learning_rate": best_params["learning_rate"],
        "weight_decay": best_params["weight_decay"],
        "batch_size": best_params["batch_size"],
    }
    print(f"Loaded study: {study.study_name}")

    # loop over latent dimensions
    for latent_dim in modes_list:
        # trains and saves a model
        #train_fieldautoencoder(latent_dim=latent_dim, model_subfolder="latent_sweep", **configuration)
        train_fieldautoencoder(latent_dim=latent_dim,model_subfolder="latent_sweep_v2",**configuration)
    
def main():
    # this is for training one specific autoencoder
    #train_fieldautoencoder(latent_dim=16, hidden_dim_1=256, hidden_dim_2=128, activation="silu",batch_size=16, learning_rate=0.000167182789476516, weight_decay=4.1083188946999205e-06, model_subfolder="testing")
    
    # this performs the hyperparameter optimization and trains one autoencoder for every entry in modes_list
    hyperparameter_optimization()
    modes_list = np.arange(4, 61, 4)
    #train_autoencoders_for_comparison(modes_list)

if __name__ == "__main__":
    main()
