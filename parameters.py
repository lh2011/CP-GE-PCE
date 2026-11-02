import numpy as np
import json

PARAMETER_RANGES = {
    "alpha":  (np.deg2rad(0), np.deg2rad(20)),
    "MH":     (1.5e-3, 12e-3),
    "MW":     (7e-3, 23e-3),
    "MAG":    (5e-3, 15e-3),
    "LW1":    (1e-3, 3e-3),
    "LW2":    (2e-3, 4e-3),
    "Theta1": (np.deg2rad(18), np.deg2rad(22)),
    "Theta2": (np.deg2rad(0.1), np.deg2rad(0.6)),
    "IApp":   (9.0, 11.0),
    "Phase0": (np.deg2rad(-5), np.deg2rad(5)),
    }

ALL_PARAMETERS = list(PARAMETER_RANGES)

# Values used when a parameter is not sampled
PARAMETER_DEFAULTS = {
    "alpha":  np.deg2rad(10),
    "MH":     7e-3,
    "MW":     19e-3,
    "MAG":    7e-3,
    "LW1":    2e-3,
    "LW2":    3e-3,
    "Theta1": np.deg2rad(21.2),
    "Theta2": np.deg2rad(0.32),
    "IApp":   10,
    "Phase0": np.deg2rad(0),
    }

def get_config_json():
    """This return a json containing the ranges and defaults so they can be read via MATLAB."""
    return json.dumps({
        "ranges": {name: [float(value) for value in bounds] for name, bounds in PARAMETER_RANGES.items()},
        "defaults": {name: float(value) for name, value in PARAMETER_DEFAULTS.items()},
    })
