%function [] = fieldPlot()
warning("This can only plot for a sampled alpha or the default one")
addpath(genpath(fullfile("packages")));

%% Get the simulation parameters
filename = "data" + filesep + "dataset_p_5_samp_623.csv";
datatable = readtable(filename);
tableNames = string(datatable.Properties.VariableNames);

%% extract the sample we want
i = 572;

%% load config file from python (parameter ranges and defaults)
repoDir = fileparts(mfilename('fullpath'));
if count(py.sys.path, repoDir) == 0
    insert(py.sys.path, int32(0), repoDir);
end
parameterModule = py.importlib.import_module('parameters');
parameterModule = py.importlib.reload(parameterModule);
config = jsondecode(char(parameterModule.get_config_json()));

% Parameter ranges: same values as in the Python sampling script
parameterRanges = config.ranges;
alpha_default = config.defaults.alpha;
IApp_default = config.defaults.IApp;
Phase0_default = config.defaults.Phase0;
IApp_min = parameterRanges.IApp(1);
IApp_max = parameterRanges.IApp(2);
phaseOffset_min = parameterRanges.Phase0(1);
phaseOffset_max = parameterRanges.Phase0(2);
alpha_min = parameterRanges.alpha(1);
alpha_max = parameterRanges.alpha(2);

%% Check which operating parameters are sampled in this dataset
sampleAlpha = ismember("alpha", tableNames);
sampleIApp = ismember("IApp", tableNames);
samplePhase0 = ismember("Phase0", tableNames);

allParameters = {'MH','MW','MAG','LW1','LW2','Theta1','Theta2'};
data_names = allParameters( ismember(allParameters, datatable.Properties.VariableNames));

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

if sampleAlpha
    RotationAngle = datatable.alpha(i);
else
    RotationAngle = alpha_default;
end
if sampleIApp
    IApp = datatable.IApp(i);
else
    IApp = IApp_default;
end
if samplePhase0
    Phase0 = datatable.Phase0(i);
else
    Phase0 = Phase0_default;
end

StatorWindings.resetExcitations();
PhaseU = EXC_CurrentSine_OPT(StatorWindings, phaseDefU, IApp, pp,  0, Phase0, "NumberWindings", NumberWindings, "PlotColor", TUDa_getColor("9c"));
PhaseV = EXC_CurrentSine_OPT(StatorWindings, phaseDefV, IApp, pp, -2/3*pi, Phase0, "NumberWindings", NumberWindings, "PlotColor", TUDa_getColor("4d"));
PhaseW = EXC_CurrentSine_OPT(StatorWindings, phaseDefW, IApp, pp,  2/3*pi, Phase0, "NumberWindings", NumberWindings, "PlotColor", TUDa_getColor("1c"));

% Coupling
Coupling = CPL_HarmonicMortaring_OPT(Rotor, 2, Stator, 2);
Coupling.setHarmonics("automatic", "SymmetryType", "antiperiodic", "DofRatio", 0.5);

Rotor.resetOptimizationResponses();
Rotor.resetOptimizationParameters();
Rotor.addOptimizationResponse("Torque");

% Phase0 in the motor contains the electrical rotor angle plus the sampled/default phase offset. Its bounds must cover the total angle.
phaseTotalMin = pp*alpha_min + phaseOffset_min;
phaseTotalMax = pp*alpha_max + phaseOffset_max;

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
% Register Stator + rotation parameter values
Stator.addOptimizationParameter("IApp", IApp, IApp_min, IApp_max);
Stator.addOptimizationParameter("Phase0", pp*RotationAngle + Phase0, phaseTotalMin, phaseTotalMax);

Coupling.addOptimizationParameter("RotationAngle", RotationAngle, alpha_min, alpha_max);

% set rotor geometry function to allow parameter variation
Rotor.setGeometryFunction(@GEO_IPMrotor_I_1);

Coupling.resetOptimizationResponses();
Coupling.addOptimizationResponse("Torque");

%% newton solver setup
Motor = OPT_System();
Motor.NewtonTolerance = 1e-15;
Motor.addElement(Rotor);
Motor.addElement(Stator);
Motor.addElement(Coupling);

Motor.resetOptimizationConstraints();
Motor.resetOptimizationObjectives();
Motor.addOptimizationObjective("Torque", 1);

%% add parameter values or default if not sampled
optNames = string({Motor.OptimizationParameters.Name});
param_now = [Motor.OptimizationParameters.Init];

% read geometric param values from sample
for iopt = 1:numel(data_names)
    name = string(data_names{iopt});
    value = datatable.(name)(i);
    param_now = setParameterByName(param_now, optNames, name, value);
end
Motor.optimizationUpdate(param_now);

%% solve motor
y0_full = zeros(Motor.NumberDOF, 1); % initial guess
[y0_full, ~] = Motor.solveStatic(0, y0_full, "Verbose", true); % full solution including bc and Jac

%% plot ground truth
plot_and_save_field(i, Rotor, Stator, "groundtruth", true);

%% PCE: load and plot
pce_data = load("data/reconstructed_fields/pce_POD_p5_order6_part4.mat");
Y_pred_pce_mat = [pce_data.Y_train_pred, pce_data.Y_val_pred, pce_data.Y_test_pred];
Y_pred_pce = Y_pred_pce_mat(:,i);
y0_full_pce = Motor.reconstructSolution(0, Y_pred_pce); 
Motor.postprocess(0, y0_full_pce);

plot_and_save_field(i, Rotor, Stator, "pce", true)

%% PCE vs. ground truth:
y0_full_diff_pce = Motor.reconstructSolution(0, y0_full(Motor.K) - Y_pred_pce); 
Motor.postprocess(0, y0_full_diff_pce);

plot_and_save_field(i, Rotor, Stator, "pce_diff", false)

%% GE-PCE: load and plot
gepce_data = load("data/reconstructed_fields/gepce_POD_p5_order6_part4.mat");
Y_pred_gepce_mat = [gepce_data.Y_train_pred, gepce_data.Y_val_pred, gepce_data.Y_test_pred];
Y_pred_gepce = Y_pred_gepce_mat(:,i);
y0_full_gepce = Motor.reconstructSolution(0, Y_pred_gepce); 
Motor.postprocess(0, y0_full_gepce);

plot_and_save_field(i, Rotor, Stator, "gepce", true)

%% PCE vs. ground truth:
y0_full_diff_gepce = Motor.reconstructSolution(0, y0_full(Motor.K) - Y_pred_gepce); 
Motor.postprocess(0, y0_full_diff_gepce);

plot_and_save_field(i, Rotor, Stator, "gepce_diff", false)

%% recompute rel errors
K_pce = pce_data.K_mat;
rel_pce = 100*sqrt(((y0_full(Motor.K) - Y_pred_pce)'*K_pce*(y0_full(Motor.K) - Y_pred_pce)) ./ (y0_full(Motor.K)'*K_pce*y0_full(Motor.K)))

K_gepce = gepce_data.K_mat;
rel_gepce = 100*sqrt(((y0_full(Motor.K) - Y_pred_gepce)'*K_gepce*(y0_full(Motor.K) - Y_pred_gepce)) ./  (y0_full(Motor.K)'*K_gepce*y0_full(Motor.K)))

%% function
function [] = plot_and_save_field(i, Rotor, Stator, suffix, plotIsolines)
    plotName = "MagneticFluxDensity";

    % Colour background: magnetic flux-density magnitude
    Rotor.plotScalar("MagneticFluxDensity", 0, "FlattenPlot", true);
    Stator.plotScalar("MagneticFluxDensity", 0, "FlattenPlot", true);

    % Use identical contour levels for rotor and stator
    if plotIsolines
            Az = [Rotor.Solution.Values; Stator.Solution.Values];
            levelLimits = [min(Az(:)), max(Az(:))];
        for region = {Rotor, Stator}
            r = region{1};
            r.plotContour("MagneticVectorPotential", 0, "FigureName", plotName, "Filled", false, "LevelLimits", levelLimits, "nLevels", 16, "LineColor", [0.65, 0.65, 0.65], "LineWidth", 0.5);
            r.plotBoundary(0, "FigureName", plotName, "ColorInterfaces", [0.3, 0.3, 0.3], "LineWidthInterfaces", 0.25);
            Az = [Rotor.Solution.Values; Stator.Solution.Values];
            levelLimits = [min(Az(:)), max(Az(:))];
        end
    end

    ax = gca;
    if plotIsolines
        clim(ax, [0, 2]);
    else
        clim(ax, [0, 0.1]);
    end
    view(ax, -60, 90);   % Rotate the displayed plot 90° counterclockwise
    axis(ax, "equal");
    axis(ax, "tight");
    axis(ax, "off");     % Hide axes, ticks, and labels
    colormap(viridis)
    cb = colorbar(ax);
    cb.Units = "normalized";
    cb.Position = [0.82, 0.10, 0.025, 0.65];
    cb.FontSize = 21;

    cb.Label.Interpreter = "latex";
    cb.Label.FontSize = 23;

    if plotIsolines
        cb.Ticks = 0:0.5:2;
        %cb.Label.String = '$|\boldsymbol{B}|\;(\mathrm{T})$';
    else
        cb.Ticks = 0:0.025:0.1;
        %cb.Label.String = '$|\boldsymbol{B}_{\mathrm{ref}}-\boldsymbol{B}_{\mathrm{pred}}|\;(\mathrm{T})$';
    end

    name_extension = sprintf("%s_%d", suffix, i);
    exportgraphics(gcf, "plots" + filesep + "Fieldplot_" + name_extension + ".png", "ContentType", "image", 'Resolution', 150)
    close all
end
