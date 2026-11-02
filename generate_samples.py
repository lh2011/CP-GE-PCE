import numpy as np
import pandas as pd
from scipy.stats import qmc
from parameters import PARAMETER_RANGES, PARAMETER_DEFAULTS
import sys
import os
current_dir = os.path.dirname(__file__) if '__file__' in locals() else os.getcwd()
target_dir = os.path.join(current_dir, '..', 'localfuns')
sys.path.insert(0, target_dir)

# choose active parameter set
#active_parameters = ["alpha","MH", "MW"]
active_parameters = ["alpha","MH", "MW", "MAG", "Theta1"]
#active_parameters = ["MH", "MW", "MAG", "LW1","LW2", "Theta1","Theta2"]
#active_parameters = ["alpha", "MH", "MW", "MAG", "LW1","LW2", "Theta1","Theta2","IApp","Phase0"]

if "alpha" not in active_parameters:
    # choose gridded rotor positions
    MECH = np.arange(0, 20, 1.5)

# sobol sampling function
def mz_sobol(partitions, ranges, names, seed):
    # partitions = 2 exponent 
    nparams = ranges.shape[0]
    
    # creates nested sobol
    sampler = qmc.Sobol(d=nparams, scramble=True, seed = seed)
    samples = sampler.random_base2(m=partitions[-1])
    samples = mz_scale(samples, ranges)
    
    # create id for points
    num_partitions = partitions.size
    part_id = np.ones((2**partitions[-1]))
    
    for ii in range(num_partitions):
        final = 2**partitions[num_partitions-(ii+1)]
        part_id[:final] = part_id[:final]*0+(num_partitions-ii-1)   
        
    # build dataframe
    df_columns = ['sobol_id'] + names
    df = pd.DataFrame(np.column_stack((part_id, samples)), columns=df_columns)
    
    # Convert 'id' to int for clarity
    df['sobol_id'] = df['sobol_id'].astype(int)
    
    return df

def mz_scale(samples, ranges):
    nparams = ranges.shape[0]
    
    for ii in range(nparams):
        minv = ranges[ii, 0]
        maxv = ranges[ii, 1]
        samples[:,ii] = samples[:,ii] * (maxv-minv)+minv  
    
    return samples

# checks if parameter combination of 3 (MW, MH, MAG) is legal
def isItLegal3(row):
    # check 1: bottom of magnet contained in rotor
    check1 = row.RRO - row.MAG -row.MH-row.RRI>0
    
    # check2: magnet not over air
    tolerance = 10 # degrees of tol for angles
    
    p4 = np.array([0.5*row.MW, row.RRO - row.MAG])
    p_outer_1 = row.RRO*np.array([np.sin(row.Theta1), np.cos(row.Theta1)])
    p_wing_1 = p_outer_1 -(p_outer_1/np.linalg.norm(p_outer_1))*row.LW1
    vector1 = p_wing_1 - p4 
    angle1 = np.rad2deg((np.atan2(vector1[1],vector1[0])+ 2*np.pi) % (2*np.pi))

    p3 = np.array([0.5*row.MW, row.RRO - row.MAG-row.MH])
    p_outer_2 = row.RRO*np.array([np.sin(row.Theta1+row.Theta2), np.cos(row.Theta1+row.Theta2)])
    p_wing_2 = p_outer_2 -(p_outer_2/np.linalg.norm(p_outer_2))*row.LW2
    vector2 = p_wing_2 - p3 
    angle2 = np.rad2deg((np.atan2(vector2[1],vector2[0])+ 2*np.pi) % (2*np.pi))

    angle3 = 90 - np.rad2deg(np.atan(p3[1]/p3[0]))
    
    check2 = (angle1<(90-tolerance))*(angle2<(90-tolerance))*(angle3<30)
    
    # check3: air not over rotor
    check3 = (row.Theta1 + row.Theta2) < np.deg2rad(30)
    
    # all legal check
    legal = check1*check2*check3

    return legal*1

def main(active_parameters):
    # Check parameter names
    unknown_parameters = set(active_parameters) - set(PARAMETER_RANGES)
    if unknown_parameters:
        raise ValueError(f"Unknown parameters: {sorted(unknown_parameters)}")
    if len(active_parameters) == 0:
        raise ValueError("At least one active parameter must be selected.")

    # Assemble names and ranges in the requested order
    names = list(active_parameters)
    ranges = np.array([PARAMETER_RANGES[name] for name in names],dtype=float,)
    dimensions = ranges.shape[0]

    # 2^m samples are part of every partition
    p_train = np.array([5,6,7])
    p_test = np.array([6])
    p_val = np.array([1])

    df_train = mz_sobol(p_train, ranges, names, 42)
    df_val = mz_sobol(p_val, ranges, names, 7)
    df_test = mz_sobol(p_test, ranges, names, 24)

    df_train.insert(0, "dataset", 0)
    df_val.insert(0, "dataset", 1)
    df_test.insert(0, "dataset", 2)

    df_all = pd.concat([df_train, df_val, df_test], ignore_index=True)

    print('dim_train =', df_train.shape)
    print('dim_val =', df_val.shape)
    print('dim_test =', df_test.shape)

    # setup legality df
    totsamples = df_all.shape[0]

    # df plus variables required for legality check
    df_general = df_all.copy()

    # Add default values for parameters that were not sampled
    for parameter, default_value in PARAMETER_DEFAULTS.items():
        if parameter not in df_general.columns:
            df_general[parameter] = default_value

    # Add fixed motor dimensions required by isItLegal3()
    df_general["RRO"] = 44e-3
    df_general["RRI"] = 16e-3  

    # %% check legality
    for ii in range(totsamples):
        # check legality status
        row = df_general.iloc[ii]
        legal = isItLegal3(row)
        df_general.at[ii, "legal"] = legal
        
    # %% alarm if illegal solution present
    if df_general["legal"].min() != 1:
        print("illegal combination found")  
        # keep only legal
        df_general = df_general[df_general["legal"] == 1]
        totsamples = df_general.shape[0]

    # Save dataset metadata and only the sampled parameters
    output_columns = ["dataset","sobol_id",*active_parameters,]

    df_all = df_general[output_columns].copy()

    df_all.insert(2, "global_id", range(1, len(df_all) + 1))

    # print how many samples
    print("train", (df_all['dataset'] == 0).sum())
    print("val", (df_all['dataset'] == 1).sum())
    print("test", (df_all['dataset'] == 2).sum())

    # Save combined dataset
    df_all.to_csv(f"data/dataset_p_{dimensions}_samp_{totsamples}.csv", index=False)
    # check how we handle alpha
    if "alpha" not in active_parameters:
        np.savetxt(f"data/MECH_ANGLE_p_{dimensions}_samp_{totsamples}.csv", MECH, delimiter=",")
    print("✅ Combined dataset saved")

def generate_n_sets(active_parameters, num_sets):
    seeds = np.array([42, 777, 2001, 123, 666, 13, 2026, 999, 1, 69])
    if num_sets > seeds.size:
        raise ValueError("Add more seeds!")

    # Check parameter names
    unknown_parameters = set(active_parameters) - set(PARAMETER_RANGES)
    if unknown_parameters:
        raise ValueError(f"Unknown parameters: {sorted(unknown_parameters)}")
    if len(active_parameters) == 0:
        raise ValueError("At least one active parameter must be selected.")

    # Assemble names and ranges in the requested order
    names = list(active_parameters)
    ranges = np.array([PARAMETER_RANGES[name] for name in names],dtype=float,)
    dimensions = ranges.shape[0]

    # 2^m samples are part of every partition
    p_train = np.array([5,6,7])
    p_test = np.array([6])
    p_val = np.array([1])

    for seed in seeds[:num_sets]:
        # sample values
        df_train = mz_sobol(p_train, ranges, names, seed)
        df_val = mz_sobol(p_val, ranges, names, seed+1)
        df_test = mz_sobol(p_test, ranges, names, seed+2)

        df_train.insert(0, "dataset", 0)
        df_val.insert(0, "dataset", 1)
        df_test.insert(0, "dataset", 2)

        df_all = pd.concat([df_train, df_val, df_test], ignore_index=True)

        print('dim_train =', df_train.shape)
        print('dim_val =', df_val.shape)
        print('dim_test =', df_test.shape)

        # setup legality df
        totsamples = df_all.shape[0]

        # df plus variables required for legality check
        df_general = df_all.copy()

        # Add default values for parameters that were not sampled
        for parameter, default_value in PARAMETER_DEFAULTS.items():
            if parameter not in df_general.columns:
                df_general[parameter] = default_value

        # Add fixed motor dimensions required by isItLegal3()
        df_general["RRO"] = 44e-3
        df_general["RRI"] = 16e-3  

        # %% check legality
        for ii in range(totsamples):
            # check legality status
            row = df_general.iloc[ii]
            legal = isItLegal3(row)
            df_general.at[ii, "legal"] = legal
            
        # %% alarm if illegal solution present
        if df_general["legal"].min() != 1:
            print("illegal combination found")  
            # keep only legal
            df_general = df_general[df_general["legal"] == 1]
            totsamples = df_general.shape[0]

        # Save dataset metadata and only the sampled parameters
        output_columns = ["dataset","sobol_id",*active_parameters,]

        df_all = df_general[output_columns].copy()

        df_all.insert(2, "global_id", range(1, len(df_all) + 1))

        # print how many samples
        print("train", (df_all['dataset'] == 0).sum())
        print("val", (df_all['dataset'] == 1).sum())
        print("test", (df_all['dataset'] == 2).sum())

        # Save combined dataset
        df_all.to_csv(f"data/cp_eval_data/dataset_p_{dimensions}_samp_{totsamples}_seed_{seed}.csv", index=False)
        # check how we handle alpha
        if "alpha" not in active_parameters:
            np.savetxt(f"data/cp_eval_data/MECH_ANGLE_p_{dimensions}_samp_{totsamples}_seed_{seed}.csv", MECH, delimiter=",")
        print(f"✅ Combined dataset saved for seed {seed}")

if __name__ == "__main__":
    #main(active_parameters)
    generate_n_sets(active_parameters, num_sets=10)
