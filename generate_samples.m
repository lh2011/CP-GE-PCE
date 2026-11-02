%% run this and change the filename to run this for one specific dataset
generate_snapshots("data/dataset_p_10_samp_2978.csv")

%% run this to generate the 10 datasets used for evaluating empirical coverage
%pattern = fullfile("data", "cp_eval_data", "dataset_p_5_samp_*_seed_*.csv");
%assert(numel(dir(pattern)) == 10, "Expected exactly 10 CP dataset files.");
%generate_snapshots(pattern, false)

function generate_snapshots(datasetPattern, savePlotData)
% Generate snapshots for one CSV or a pattern matching several CSVs.
if nargin < 1
    datasetPattern = fullfile("data", "dataset_p_5_samp_623.csv");
end
if nargin < 2
    savePlotData = false;
end

addpath(genpath(fullfile("packages", "amiga-fork")));
files = dir(datasetPattern);
assert(~isempty(files), 'No dataset files match "%s".', datasetPattern);

% load config file from python (parameter ranges and defaults)
repoDir = fileparts(mfilename('fullpath'));
if count(py.sys.path, repoDir) == 0
    insert(py.sys.path, int32(0), repoDir);
end
parameterModule = py.importlib.import_module('parameters');
parameterModule = py.importlib.reload(parameterModule);
config = jsondecode(char(parameterModule.get_config_json()));

[~, order] = sort({files.name});
files = files(order);
for k = 1:numel(files)
    dataname = fullfile(files(k).folder, files(k).name);
    fprintf("\nDataset %d/%d: %s\n", k, numel(files), dataname);

    simulateOneDataset(dataname, savePlotData, config);
end
end

function simulateOneDataset(dataname, savePlotData, config)
datatable = readtable(dataname);
tableNames = string(datatable.Properties.VariableNames);

[fileDir, fileStem] = fileparts(dataname);
fileStem = string(fileStem);
assert(startsWith(fileStem, "dataset_"));

sampleAlpha = ismember("alpha", tableNames);
if sampleAlpha
    alpha_samples = datatable.alpha;
    num_angles = 1;
    angle_mode = "sampled";
else
    mechname = fullfile(fileDir, replace(fileStem, "dataset_", "MECH_ANGLE_") + ".csv");
    assert(isfile(mechname), "Missing angle file: %s", mechname);
    mechtable = readmatrix(mechname);
    mechtable = mechtable(:);
    num_angles = numel(mechtable);
    angle_mode = "legacy_grid";
end

num_samps = height(datatable);

% Check which operating parameters are sampled in this dataset
sampleIApp = ismember("IApp", tableNames);
samplePhase0 = ismember("Phase0", tableNames);
fprintf("Sample IApp: %d\n", sampleIApp);
fprintf("Sample Phase0: %d\n", samplePhase0);

allParameters = {'MH','MW','MAG','LW1','LW2','Theta1','Theta2'};
data_names = allParameters( ismember(allParameters, datatable.Properties.VariableNames));

%% Parameter ranges: same values as in the Python sampling script
parameterRanges = config.ranges;
IApp_default = config.defaults.IApp;
Phase0_default = config.defaults.Phase0;
IApp_min = parameterRanges.IApp(1);
IApp_max = parameterRanges.IApp(2);
phaseOffset_min = parameterRanges.Phase0(1);
phaseOffset_max = parameterRanges.Phase0(2);
alpha_min = parameterRanges.alpha(1);
alpha_max = parameterRanges.alpha(2);

%% setup motor
% Rotor
% Geometry
[srfRotor, regionsRotor, boundariesRotor, ~, mag] = GEO_IPMrotor_I_1(preview=false);
% Discretization
Rotor = IGA_RegionMagnetic_OPT([5, 5], [2,2], 0.1);
Rotor.importSurface(srfRotor, regionsRotor, boundariesRotor);
% Materials
Rotor.resetMaterials();
RotorCore = MAT_M27(Rotor, 1);
RotorAir = MAT_Air(Rotor, 2);
RotorMagnet = MAT_Magnet_OPT(Rotor, 3, mag, "Remanence", 1.4);
% Boundary Conditions
Rotor.resetBoundaryConditions()
RotorBCdirichlet = IGA_BCdirichlet(Rotor, 1);
RotorBCap = IGA_BCantiPeriodic(Rotor, 3, 4);

% Stator
% Geometry
[srfStator, regionsStator, boundariesStator, ~] = GEO_DISstator_1("preview", false, ...
    "Symmetry", 6, "SLOTS", 36, ...
    "HSO", 0.64e-3, "HY", 8.25e-3, ...
    "RAG", 44.5e-3, "SRI", 45e-3, "SRO", 67.5e-3, ...
    "SR1", 1.2e-3, "WSO", 2.3e-3, "WT", 4e-3, "LumpWindings", false);
% Discretization
Stator = IGA_RegionMagnetic_OPT([3, 3], [2,2], 0.1);
Stator.importSurface(srfStator, regionsStator, boundariesStator);
% Materials
StatorCore = MAT_M27(Stator, 1);
StatorAir = MAT_Air(Stator, 2);
StatorWindings = MAT_Copper(Stator, 3:14);
% Boundary Conditions
Stator.resetBoundaryConditions()
IGA_BCdirichlet(Stator,1);
IGA_BCantiPeriodic(Stator, 3, 4);

% Excitation
phaseDefU(1).Regions = 3;       phaseDefU(1).Direction = 1;
phaseDefU(2).Regions = 4;       phaseDefU(2).Direction = 1;
phaseDefU(3).Regions = 6;       phaseDefU(3).Direction = 1;
phaseDefU(4).Regions = 13;      phaseDefU(4).Direction = -1;

phaseDefV(1).Regions = 9;       phaseDefV(1).Direction = 1;
phaseDefV(2).Regions = 11;      phaseDefV(2).Direction = 1;
phaseDefV(3).Regions = 12;      phaseDefV(3).Direction = 1;
phaseDefV(4).Regions = 14;      phaseDefV(4).Direction = 1;

phaseDefW(1).Regions = 5;       phaseDefW(1).Direction = -1;
phaseDefW(2).Regions = 7;       phaseDefW(2).Direction = -1;
phaseDefW(3).Regions = 8;       phaseDefW(3).Direction = -1;
phaseDefW(4).Regions = 10;      phaseDefW(4).Direction = -1;

NumberWindings = 12;
pp = 3;

StatorWindings.resetExcitations();
PhaseU = EXC_CurrentSine_OPT(StatorWindings, phaseDefU, IApp_default, pp,  0, Phase0_default, "NumberWindings", NumberWindings, "PlotColor", TUDa_getColor("9c"));
PhaseV = EXC_CurrentSine_OPT(StatorWindings, phaseDefV, IApp_default, pp, -2/3*pi, Phase0_default, "NumberWindings", NumberWindings, "PlotColor", TUDa_getColor("4d"));
PhaseW = EXC_CurrentSine_OPT(StatorWindings, phaseDefW, IApp_default, pp,  2/3*pi, Phase0_default, "NumberWindings", NumberWindings, "PlotColor", TUDa_getColor("1c"));

% Coupling
Coupling = CPL_HarmonicMortaring_OPT(Rotor, 2, Stator, 2);
Coupling.setHarmonics("automatic", "SymmetryType", "antiperiodic", "DofRatio", 0.5);

Rotor.resetOptimizationResponses();
Rotor.resetOptimizationParameters();
Rotor.addOptimizationResponse("Torque");

%% check admissibility of parameters and create matrices
currentAmplitudes = repmat(IApp_default, num_samps, 1);
phaseOffsets = repmat(Phase0_default, num_samps, 1);

if sampleIApp
    currentAmplitudes = datatable.IApp;
end
if samplePhase0
    phaseOffsets = datatable.Phase0;
end

% check that all samples angles lie within the boundary
if sampleAlpha
    assert(all(alpha_samples >= alpha_min & alpha_samples <= alpha_max), 'Some sampled alpha values lie outside alpha_min/max.');
end

% Phase0 in the motor contains the electrical rotor angle plus the sampled/default phase offset. Its bounds must cover the total angle.
phaseTotalMin = pp*alpha_min + phaseOffset_min;
phaseTotalMax = pp*alpha_max + phaseOffset_max;

if sampleAlpha
    initialRotationAngle = alpha_samples(1);
else
    initialRotationAngle = mechtable(1);
end

%% Add each active parameter as an optimization parameter
for iopt = 1:numel(data_names)
    name = string(data_names{iopt});
    val = datatable.(name)(1);
    range = parameterRanges.(name);
    valm  = range(1);
    valM  = range(2);
    Rotor.addOptimizationParameter(name, val, valm, valM);
end

% Optimization responses
Stator.resetOptimizationResponses();
Stator.resetOptimizationParameters();
Stator.addOptimizationResponse("Torque");
% Always register the operating parameters. If they are not sampled, their default values are assigned in the solve loop.
Stator.addOptimizationParameter("IApp", IApp_default, IApp_min, IApp_max);
Stator.addOptimizationParameter("Phase0", pp*initialRotationAngle + Phase0_default, phaseTotalMin, phaseTotalMax);

Coupling.addOptimizationParameter("RotationAngle", initialRotationAngle, alpha_min, alpha_max);

% set rotor geometry function to allow parameter variation
Rotor.setGeometryFunction(@GEO_IPMrotor_I_1);

Coupling.resetOptimizationResponses();
Coupling.addOptimizationResponse("Torque");

% newton solver setup
Motor = OPT_System();
Motor.NewtonTolerance = 1e-15;
Motor.addElement(Rotor);
Motor.addElement(Stator);
Motor.addElement(Coupling);

Motor.resetOptimizationConstraints();
Motor.resetOptimizationObjectives();
Motor.addOptimizationObjective("Torque", 1);

%% prepare variables for solve cycle
% first guess = 0
y0_full = zeros(Motor.NumberDOF, 1); % initial guess
[y0_full, ~] = Motor.solveStatic(0, y0_full, "Verbose", true); % full solution including bc and Jac
DOFs_num = numel(Motor.K);

%% preallocate
state_store = zeros(num_angles, DOFs_num, num_samps);
state_der_store = zeros(num_angles, DOFs_num, num_samps, Motor.NumberOptimizationParameters);

global_store = zeros(num_angles, num_samps);
global_der_store = zeros(num_angles, num_samps,  Motor.NumberOptimizationParameters);

%% solve cycle
solveTimer = tic;
% create parameter vector
optNames = string({Motor.OptimizationParameters.Name});
derivative_names = cellstr(optNames);

idxRotation = find(optNames == "RotationAngle", 1);
idxPhase0   = find(optNames == "Phase0", 1);

if sampleAlpha
    % The transformed RotationAngle column represents the derivative with respect to the PCE input alpha.
    derivative_names{idxRotation} = 'alpha';
end

param_now = [Motor.OptimizationParameters.Init];
disp(table((1:numel(optNames))', optNames', 'VariableNames', {'Index', 'Parameter'}));

rotation_angle_store = zeros(num_angles, num_samps);

times_pce = zeros(num_samps, 1);
times_gepce = zeros(num_samps, 1);

% loop over rotation angles
for jj = 1:num_angles
    % loop over samples
    for ii = 1:num_samps
        tic();
        % print iteration 
        disp(num2str(jj) + " / " + num2str(num_angles) + " (" + num2str(ii) + " / " + num2str(num_samps) + ")")
        
        % we need to distinguish here to allow for the old case where we use a table of angles for every parameter combination and the new case where the angle is sampled aswell
        if sampleAlpha
            rotationAngle = alpha_samples(ii);
        else
            rotationAngle = mechtable(jj);
        end
        rotation_angle_store(jj, ii) = rotationAngle;

        % read geometric param values from sample
        for iopt = 1:numel(data_names)
            name = string(data_names{iopt});
            value = datatable.(name)(ii);
            param_now = setParameterByName(param_now, optNames, name, value);
        end
        param_now = setParameterByName(param_now, optNames, "IApp", currentAmplitudes(ii));
        totalElectricalPhase = pp * rotationAngle + phaseOffsets(ii);
        param_now = setParameterByName(param_now, optNames, "Phase0", totalElectricalPhase);
        param_now = setParameterByName(param_now, optNames,"RotationAngle", rotationAngle);
        % update params
        Motor.optimizationUpdate(param_now);
        
        % field value
        [y0_full, Jac] = Motor.solveStatic(0, y0_full, "Verbose", false); % full solution including bc and Jac
        y0 = y0_full(Motor.K);

        times_pce(ii) = toc();
    
        % field derivative
        dedx = Motor.partial_e_partial_x(0, param_now, y0_full); %full system
        dedxBC = dedx(Motor.K, :) + Motor.G'*(dedx(Motor.R, :));
        % dydx reduced with boundary conditions
        dydxBC = -Jac\dedxBC; % this is field derivative

        % Application of the chain rule because we use 
        % phase offset = pp * alpha + phaseOffset_0
        if sampleAlpha
            dydxBC(:, idxRotation) = dydxBC(:, idxRotation) + pp*dydxBC(:, idxPhase0);
        end

        times_gepce(ii) = toc();

        % torque and derivative computation
        [T, Tder] = Motor.response_function(param_now);
        if sampleAlpha
            Tder(idxRotation) = Tder(idxRotation) + pp*Tder(idxPhase0);
        end
    
        % store results
        global_store(jj, ii) = T;
        global_der_store(jj,ii,:) = Tder;
        state_store(jj,:, ii) = y0;
        state_der_store(jj,:, ii, :) = dydxBC;
    end
end
gentime = toc(solveTimer);

disp(mean(times_pce(2:11)))
disp(mean(times_gepce(2:11)))

%% save solutions
filename = fullfile(fileDir, fileStem + ".mat");
% save solution
save(filename,"global_store", "global_der_store","state_store", "state_der_store", "rotation_angle_store", "angle_mode", "derivative_names", "gentime", "-v7.3")
% save stuff 4 plots
if savePlotData
    plotFilename = fullfile(fileDir, replace(fileStem, "dataset_", "plot_data_") + ".mat");
    save(plotFilename, "Motor", "Rotor", "Stator", "Coupling");
end

%% Compute error matrices at the mean parameter combination
% Start with the current/default optimization parameter values
p_mean = [Motor.OptimizationParameters(:).Init];
optNames = string({Motor.OptimizationParameters.Name});

% Set every active geometric parameter to its sample mean
for iopt = 1:numel(data_names)
    name = string(data_names{iopt});
    idx = optNames == name;
    assert(any(idx),'Parameter "%s" was not found in Motor.OptimizationParameters.', name);
    p_mean(idx) = mean(datatable.(name));
end

% Set operating parameters used for the scaling matrices
if sampleAlpha
    rotationForMatrices = mean(datatable.alpha);  % alpha stored in radians
    nAngleBlocks = 1;
else
    % Keep exactly for reproducing the previous results
    rotationForMatrices = deg2rad(10);
    nAngleBlocks = numel(mechtable);
end
idx = optNames == "RotationAngle";
if any(idx)
    p_mean(idx) = rotationForMatrices;
end

p_mean = setParameterByName(p_mean, optNames, "IApp", mean(currentAmplitudes));
phaseForMatrices = pp * rotationForMatrices + mean(phaseOffsets);
p_mean = setParameterByName(p_mean, optNames, "Phase0", phaseForMatrices);

disp("Mean parameter combination:")
for ii = 1:numel(optNames)
    fprintf("%s: %.8g\n", optNames(ii), p_mean(ii));
end

% Update geometry and system matrices
Motor.optimizationUpdate(p_mean);
% Solve once so that the system and coupling matrices correspond to p_mean, this is required due to the inner workings of amiga to update internal stuff
y0_full = zeros(Motor.NumberDOF, 1);
Motor.solveStatic(0, y0_full, "Verbose", true);

%% Assemble mass and stiffness matrices
% rotor mesh and space
RotorMesh = MESH_2D(Rotor.MeshesEval);
RotorSpace = SPACE_2D(Rotor.SpacesEval, Rotor.Spaces.gnum, 1);
% stator mesh and space
StatorMesh = MESH_2D(Stator.MeshesEval);
StatorSpace = SPACE_2D(Stator.SpacesEval, Stator.Spaces.gnum, 1);

% extract indices
rotorIndices = Motor.Elements(1).Indices;
statorIndices = Motor.Elements(2).Indices;

% Full stiffness matrix
K = sparse(Motor.NumberDOF, Motor.NumberDOF);
K(rotorIndices, rotorIndices) = int_gradu_gradv(RotorSpace, RotorSpace, RotorMesh);
K(statorIndices, statorIndices) = int_gradu_gradv(StatorSpace, StatorSpace, StatorMesh);

% Full mass matrix
M = sparse(Motor.NumberDOF, Motor.NumberDOF);
M(rotorIndices, rotorIndices) = int_u_v(RotorSpace, RotorSpace, RotorMesh);
M(statorIndices, statorIndices) = int_u_v(StatorSpace, StatorSpace, StatorMesh);

%% Reduce matrices according to boundary and coupling constraints
k = Motor.K;
r = Motor.R;
G = Motor.G;
Mred = M(k,k) + M(k,r) * G + G' * M(r,k) + G' * M(r,r) * G;
Kred = K(k,k) + K(k,r) * G + G' * K(r,k) + G' * K(r,r) * G;

%% Create block matrix for all rotor positions
if nAngleBlocks == 1
    Kred_big = Kred;
else
    Kcells = repmat({Kred}, 1, nAngleBlocks);
    Kred_big = blkdiag(Kcells{:});
end

%% Save error matrices
matrixFilename = fullfile(fileDir, replace(fileStem, "dataset_", "error_matrices_") + ".mat");
save(matrixFilename, "Mred", "Kred", "Kred_big", "p_mean", "optNames", "-v7");
fprintf("Saved error matrices to %s\n", matrixFilename);
end

function parameters = setParameterByName(parameters, parameterNames, name, value)
    % find index for parameter
    idx = find(parameterNames == string(name));
    assert(~isempty(idx), 'Optimization parameter "%s" was not found.', name);
    assert(isscalar(idx), 'Optimization parameter "%s" occurs more than once.', name);
    % assign value
    parameters(idx) = value;
end
