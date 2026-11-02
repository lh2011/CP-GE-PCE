classdef  OPT_System < EQS_System & OPT_Element
    properties
        OptimizationObjectives, NumberOptimizationObjectives;
        OptimizationConstraints, NumberOptimizationConstraints;
        LowerBounds, UpperBounds;
        lastx;      % last design state of updated system
        lasty;      % last solution as starting values for nonlinear system
        StartTime=tic();
        xsave1, rsave1; % (Re)use responses for both objective and constraints
        xsave2, rsave2, gradrsave2;   
        OptimizationHistory;
    end

    methods
        function obj = OPT_System()
            obj.resetOptimizationObjectives();
            obj.resetOptimizationParameters();
            obj.resetOptimizationConstraints();
        end

        function addOptimizationObjective(obj, name, scaling)
            if ~exist("scaling", "var")
                scaling = 1;
            end
            pos = numel(obj.OptimizationObjectives) + 1;
            if pos > 1
                assert(~any(name == vertcat(obj.OptimizationObjectives.Name)), "This objective has already been added to the system!");
            end
            obj.OptimizationObjectives(pos).Name = name;
            obj.OptimizationObjectives(pos).Scaling = scaling;
            obj.NumberOptimizationObjectives = obj.NumberOptimizationObjectives + 1;
            % Go through elements, check responses, add Indices to Responses
            obj.updateOptimizationResponses();
            obj.updateOptimizationParameters();
        end

        function addOptimizationConstraint(obj, name, type, refValue, scaling)
            assert(any(type == ["<", ">", "="]));
            if ~exist("scaling", "var")
                scaling = 1;
            end
            pos = numel(obj.OptimizationConstraints) + 1;
            if pos > 1
                assert(~any(name == vertcat(obj.OptimizationConstraints.Name)), "This constraint has already been added to the system!");
            end
            obj.OptimizationConstraints(pos).Name = name;
            obj.OptimizationConstraints(pos).RefValue = refValue;
            obj.OptimizationConstraints(pos).Type = type;
            obj.OptimizationConstraints(pos).Scaling = scaling;
            obj.NumberOptimizationConstraints = obj.NumberOptimizationConstraints + 1;
            % Go through elements, check responses, add Indices to Responses
            obj.updateOptimizationResponses();
            obj.updateOptimizationParameters();
        end

        function x_opt = optimize(obj, x0, newOptions)
            if ~exist("x0", "var") || isempty(x0)
                x0norm = (vertcat(obj.OptimizationParameters.Init) - vertcat(obj.OptimizationParameters.Min))./(vertcat(obj.OptimizationParameters.Max) - vertcat(obj.OptimizationParameters.Min));
            else
                x0norm = (x0 - vertcat(obj.OptimizationParameters.Min))./(vertcat(obj.OptimizationParameters.Max) - vertcat(obj.OptimizationParameters.Min));
            end

            obj.OptimizationHistory = struct('iteration', {}, 'funccount', {}, 'fval', {}, 'Responses', {}, 'Objectives', {}, 'Constraints', {}, 'Parameters', {}, 'TimeElapsed', {});
            StartOptTime = string(datetime('now', 'Format', "yyyy-MM-dd-HH-mm"));
            obj.SaveDirectory = "plots/"+StartOptTime;
            mkdir(obj.SaveDirectory);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    obj.Elements(iEl).Element.SaveDirectory = obj.SaveDirectory;
                end
            end
            
            problem.solver = 'fmincon';
            problem.objective = @(x) obj.objective_function(x);
            problem.nonlcon = @(x) obj.constraint_function(x);
            problem.x0 = x0norm;
            problem.Aeq = [];
            problem.beq = [];
            problem.ub = ones(obj.NumberOptimizationParameters, 1);
            problem.lb = zeros(obj.NumberOptimizationParameters, 1);
            % Standard Options, can be overwritten
            standardOptions = optimoptions(@fmincon, ... 
                'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
                'Display', 'iter-detailed', ...
                'Algorithm', 'interior-point', ...
                'MaxIterations', 500, ...
                'OutputFcn', @(a, b, c) obj.outputFunction(a, b, c), ...
                'PlotFcn', @(a, b, c) obj.plotFunction(a, b, c), ... % 'PlotFcn', @optimplotfval, ...
                'StepTolerance', 1e-8, ... 
                'OptimalityTolerance', 1e-6, ...
                'ConstraintTolerance', 1e-6, ...                
                'ScaleProblem', false, ... % scale manually (otherwise only normalized by initial value)
                'HonorBounds', false);   % do not reset design variables that are at the bounds

            % Overwrite or extend the standard options
            if exist("newOptions", "var")
                field_names = fieldnames(newOptions);
                for iNewOpts = 1:numel (field_names)
                    standardOptions.(field_names{iNewOpts}) = newOptions.(field_names{iNewOpts});
                end
            end
            problem.options = standardOptions;
            
            % start optimization
            obj.StartTime = tic();
            x_opt_norm = fmincon(problem);

            % Rescale to physical values
            x_opt = x_opt_norm.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;
        end

        % Objective function for optimization
        function [f, gradf] = objective_function(obj, x_norm)
            % Rescale values from 0-1 to real values
            x = x_norm.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            if nargout == 2
                % check if responses and derivatives have been evaluated recently
                if ~isempty(obj.xsave2) && all(obj.xsave2 == x)
                    r = obj.rsave2;
                    gradr = obj.gradrsave2;
                else
                    [r, gradr] = obj.response_function(x);
                    obj.xsave2 = x;
                    obj.rsave2 = r;
                    obj.gradrsave2 = gradr;
                end
            else
                % check if responses have been evaluated recently
                if ~isempty(obj.xsave1) && all(obj.xsave1 == x)
                    r = obj.rsave1;
                elseif ~isempty(obj.xsave2) && all(obj.xsave2 == x)
                    r = obj.rsave2;
                else
                    r = obj.response_function(x);
                    obj.xsave1 = x;
                    obj.rsave1 = r;
                end
            end
            % see which responses are objectives
            f = 0;
            gradf = zeros(obj.NumberOptimizationParameters, 1);
            for iObjective = 1:numel(obj.OptimizationObjectives)
                ind_obj = obj.OptimizationObjectives(iObjective).Index; 
                scale = obj.OptimizationObjectives(iObjective).Scaling; 
                f = f + r(ind_obj)*scale;

                if nargout == 2
                    gradf = gradf + gradr(ind_obj, :)'*scale;
                end
            end
            % Rescale gradient of f because of parameter normalization
            gradf = gradf .* (obj.UpperBounds-obj.LowerBounds);
        end

        function [c, ceq, gradc, gradceq] = constraint_function(obj, x_norm)
            if numel(obj.OptimizationConstraints) == 0
                c = []; ceq = []; gradc = []; gradceq = [];
                return
            end
            % Rescale values from 0-1 to real values
            x = x_norm.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            if nargout > 2
                % check if responses and derivatives have been evaluated recently
                if ~isempty(obj.xsave2) && all(obj.xsave2 == x)
                    r = obj.rsave2;
                    gradr = obj.gradrsave2;
                else
                    [r, gradr] = obj.response_function(x);
                    obj.xsave2 = x;
                    obj.rsave2 = r;
                    obj.gradrsave2 = gradr;
                end
            else
                % check if responses have been evaluated recently
                if ~isempty(obj.xsave1) && all(obj.xsave1 == x)
                    r = obj.rsave1;
                elseif ~isempty(obj.xsave2) && all(obj.xsave2 == x)
                    r = obj.rsave2;
                else
                    r = obj.response_function(x);
                    obj.xsave1 = x;
                    obj.rsave1 = r;
                end
            end

            % Make values and derivatives for all constraints, do not
            % distinguish inequality and equality
            allc = zeros(obj.NumberOptimizationConstraints, 1);
            allgradc = zeros(obj.NumberOptimizationParameters, obj.NumberOptimizationConstraints);

            for iConstraint = 1:numel(obj.OptimizationConstraints)
                ind_constr = obj.OptimizationConstraints(iConstraint).Index; 
                scale = obj.OptimizationConstraints(iConstraint).Scaling;
                value = obj.OptimizationConstraints(iConstraint).RefValue;
                type = obj.OptimizationConstraints(iConstraint).Type;

                if type == "<" || type == "="
                    allc(iConstraint) = (r(ind_constr) - value)*scale;
                    if nargout > 2
                        allgradc(:, iConstraint) = gradr(ind_constr, :)*scale;
                    end
                else
                    allc(iConstraint) = -(r(ind_constr) - value)*scale;
                    if nargout > 2
                        allgradc(:, iConstraint) = -gradr(ind_constr, :)*scale;
                    end
                end
            end
            % Rescale gradient of c because of parameter normalization
            allgradc = allgradc .* (obj.UpperBounds-obj.LowerBounds);
            % get right indices for equality and inequality constraints
            ind_eq = find(vertcat(obj.OptimizationConstraints.Type) == "=");
            ind_ineq = union(find(vertcat(obj.OptimizationConstraints.Type) == ">"), ...
                find(vertcat(obj.OptimizationConstraints.Type) == "<"));
            c = allc(ind_ineq);
            ceq = allc(ind_eq);
            if nargout > 2
                gradc = allgradc(:, ind_ineq);
                gradceq = allgradc(:, ind_eq);
            end
        end
        
        function [r, gradr] = response_function(obj, x, t)
            if ~exist("t", "var")
                t = 0;
            end
            mytimer("FirstOrder", "start");
            mytimer("FirstOrder-Build", "start");
            % try
                obj.optimizationUpdate(x);
            % catch
                % r = NaN;
                % gradr = Nan;
                % return
            % end
            mytimer("FirstOrder-Build", "stop");
            if isempty(obj.lasty)
                y0 = zeros(obj.NumberDOF, 1);
            else
                y0 = obj.lasty;
            end
            mytimer("FirstOrder-Newton", "start");
            [y, Jac] = obj.solveStatic(t, y0, "Verbose", false);
            mytimer("FirstOrder-Newton", "stop");
            obj.lasty = y;
            mytimer("FirstOrder-Responses", "start");
            r = obj.responses(t, x, y);
            mytimer("FirstOrder-Responses", "stop");

            if nargout == 2
                % % adjoint equation
                % mytimer("FirstOrder-Adjoint", "start");
                % drdx = obj.partial_r_partial_x(t, x, y);
                % drdy = obj.partial_r_partial_y(t, x, y);
                % drdyBC = drdy(:, obj.K) + drdy(:, obj.R)*obj.G;
                % zk = -(Jac'\drdyBC');
                % % reconstruct full adjoint
                % zr = obj.G*zk;
                % z = zeros(obj.NumberDOF, obj.NumberOptimizationResponses);
                % z(obj.K, :) = zk;
                % z(obj.R, :) = zr;
                % mytimer("FirstOrder-Adjoint", "stop");
                % mytimer("FirstOrder-Lx", "start");
                % Lx = obj.function_Lx(t, x, y, z);
                % gradr = drdx + Lx;
                % mytimer("FirstOrder-Lx", "stop");

                mytimer("FirstOrder-Gradient", "start");
                mytimer("Building-Gradients", "start");
                drdx = obj.partial_r_partial_x(t, x, y);
                dedx = obj.partial_e_partial_x(t, x, y);
                drdy = obj.partial_r_partial_y(t, x, y);
                mytimer("Building-Gradients", "stop");
                % Applying boundary conditions
                dedxBC = dedx(obj.K, :) + obj.G'*(dedx(obj.R, :));
                drdyBC = drdy(:, obj.K) + drdy(:, obj.R)*obj.G;
                mytimer("FirstOrder-Gradient", "stop");

                mytimer("FirstOrder-GradientSolve", "start");
                if obj.NumberOptimizationObjectives > obj.NumberOptimizationParameters
                    gradr = drdx - drdyBC*(Jac\dedxBC);        % direct method
                else
                    gradr = drdx - (Jac'\drdyBC')'*dedxBC;     % adjoint method
                end
                mytimer("FirstOrder-GradientSolve", "stop");
            end
            obj.lastx = x;
            mytimer("FirstOrder", "stop");
        end

        function updateOptimizationParameters(obj)
            % if system: run again systems that elements get updated first
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_System")
                    obj.Elements(iEl).Element.updateOptimizationParameters();
                end
            end
            % Find all defined parameters for system and added elements
            parameters = vertcat(obj.OptimizationParameters.Name);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    elParameters = vertcat(obj.Elements(iEl).Element.OptimizationParameters.Name);
                    parameters = unique([parameters; elParameters], 'stable');
                end
            end
            % Find indices for parameters of added elements
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    elParameters = vertcat(obj.Elements(iEl).Element.OptimizationParameters.Name);
                    if isempty(elParameters) 
                        indices = [];
                    else
                        [~, ~, indices] = intersect(elParameters, parameters, 'stable');
                    end
                    obj.Elements(iEl).IndicesParameters = reshape(indices, 1, []);
                end
            end
            % Update System Parameters
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    for iIndex = 1:numel(obj.Elements(iEl).IndicesParameters)
                        index = obj.Elements(iEl).IndicesParameters(iIndex);
                        obj.OptimizationParameters(index).Name = obj.Elements(iEl).Element.OptimizationParameters(iIndex).Name;
                        obj.OptimizationParameters(index).Init = obj.Elements(iEl).Element.OptimizationParameters(iIndex).Init;
                        obj.OptimizationParameters(index).Min = obj.Elements(iEl).Element.OptimizationParameters(iIndex).Min;
                        obj.OptimizationParameters(index).Max = obj.Elements(iEl).Element.OptimizationParameters(iIndex).Max;
                    end
                end
            end
            obj.NumberOptimizationParameters = numel(obj.OptimizationParameters);
            obj.LowerBounds = vertcat(obj.OptimizationParameters.Min);
            obj.UpperBounds = vertcat(obj.OptimizationParameters.Max);
        end

        function updateOptimizationResponses(obj)
            % if system: run again systems such that elements get updated first
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_System")
                    obj.Elements(iEl).Element.updateOptimizationResponses();
                end
            end
            % Find all defined responses first
            responses = vertcat(obj.OptimizationResponses.Name);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    elResponses= vertcat(obj.Elements(iEl).Element.OptimizationResponses.Name);
                    responses = unique([responses; elResponses], 'stable');
                end
            end
            % Find indices for all responses
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    elResponses = vertcat(obj.Elements(iEl).Element.OptimizationResponses.Name);
                    if isempty(elResponses)
                        indices = [];
                    else
                        [~, ~, indices] = intersect(elResponses, responses, 'stable');
                    end
                    obj.Elements(iEl).IndicesResponses = reshape(indices, 1, []);
                end
            end
            % Update system responses
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    for iIndex = 1:numel(obj.Elements(iEl).IndicesResponses)
                        index = obj.Elements(iEl).IndicesResponses(iIndex);
                        obj.OptimizationResponses(index).Name = obj.Elements(iEl).Element.OptimizationResponses(iIndex).Name;
                    end
                end
            end
            obj.NumberOptimizationResponses = numel(obj.OptimizationResponses);
            % Update system objectives and constraints
            for iObjectives = 1:numel(obj.OptimizationObjectives)
                name = obj.OptimizationObjectives(iObjectives).Name;
                responseNames = vertcat(obj.OptimizationResponses.Name);
                obj.OptimizationObjectives(iObjectives).Index = find(name == responseNames);
                if isempty(obj.OptimizationObjectives(iObjectives).Index)
                    error('You specified the objective %s, which has not been set as response before', name)
                end
            end
            for iConstraint = 1:numel(obj.OptimizationConstraints)
                name = obj.OptimizationConstraints(iConstraint).Name;
                responseNames = vertcat(obj.OptimizationResponses.Name);
                obj.OptimizationConstraints(iConstraint).Index = find(name == responseNames);
                if isempty(obj.OptimizationConstraints(iConstraint).Index)
                    error('You specified the constraint %s, which has not been set as response before', name)
                end
            end
        end

        % Update all elements in the optimization 
        function optimizationUpdate(obj, x)
            if isempty(obj.lastx) || ~all(x == obj.lastx)
                % reset precalculated linear matrices for new designs
                obj.Klin = [];
                obj.Flin = [];
            end
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    obj.Elements(iEl).Element.optimizationUpdate(x(ind_x));
                end
            end
        end

        function r = responses(obj, t, x, y)
            r = zeros(obj.NumberOptimizationResponses, 1);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    ind_r = obj.Elements(iEl).IndicesResponses;
                    % ind_y = obj.Elements(iEl).Indices;
                    r_iel = obj.Elements(iEl).Element.responses(t, x_iel, y_iel);
                    r(ind_r) = r(ind_r) + r_iel;
                end
            end
        end

        function Lx = function_Lx(obj, t, x, y, z)
            Lx = zeros(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_y = obj.Elements(iEl).Indices;
                    y_iel = y(ind_y);
                    z_iel = z(ind_y, :);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    ind_r = obj.Elements(iEl).IndicesResponses;
                    
                    Lx_iel = obj.Elements(iEl).Element.function_Lx(t, x_iel, y_iel, z_iel);
                    Lx(ind_r, ind_x) = Lx(ind_r, ind_x) + Lx_iel;
                end
            end
        end

        function drdy = partial_r_partial_y(obj, t, x, y)
            drdy = zeros(obj.NumberOptimizationResponses, obj.NumberDOF);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_y = obj.Elements(iEl).Indices;
                    y_iel = y(ind_y);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    ind_r = obj.Elements(iEl).IndicesResponses;
                    
                    drdy_iel = obj.Elements(iEl).Element.partial_r_partial_y(t, x_iel, y_iel);
                    drdy(ind_r, ind_y) = drdy(ind_r, ind_y) + drdy_iel;
                end
            end
        end

        function drdx = partial_r_partial_x(obj, t, x, y)
            drdx = zeros(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    ind_r = obj.Elements(iEl).IndicesResponses;
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    drdx_iel = obj.Elements(iEl).Element.partial_r_partial_x(t, x_iel, y_iel);
                    drdx(ind_r, ind_x) = drdx(ind_r, ind_x) + drdx_iel;
                end
            end
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            dedx = zeros(obj.NumberDOF, obj.NumberOptimizationParameters);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    dedx_iel = obj.Elements(iEl).Element.partial_e_partial_x(t, x_iel, y_iel);
                    dedx(ind_iel, ind_x) = dedx(ind_iel, ind_x) + dedx_iel;
                end
            end
        end

        function stop = outputFunction(obj, xnorm, optimValues, state)
            stop = false;
            x = xnorm .* (obj.UpperBounds - obj.LowerBounds) + obj.LowerBounds;
            if ~isempty(obj.xsave1) && all(obj.xsave1 == x)
                r = obj.rsave1;
            elseif ~isempty(obj.xsave2) && all(obj.xsave2 == x)
                r = obj.rsave2;
            else
                r = obj.response_function(x);
                obj.xsave1 = x;
                obj.rsave1 = r;
            end
            % Store responses in Output struct
            pos = numel(obj.OptimizationHistory) + 1;
            obj.OptimizationHistory(pos).funccount = optimValues.funccount;
            obj.OptimizationHistory(pos).fval = optimValues.fval;
            obj.OptimizationHistory(pos).iteration = optimValues.iteration;
            obj.OptimizationHistory(pos).stepsize = optimValues.stepsize;
            obj.OptimizationHistory(pos).firstorderopt = optimValues.firstorderopt;
            obj.OptimizationHistory(pos).constrviolation = optimValues.constrviolation;
            obj.OptimizationHistory(pos).trustregionradius = optimValues.trustregionradius;
            obj.OptimizationHistory(pos).Responses = obj.OptimizationResponses;
            for iResponse = 1:obj.NumberOptimizationResponses
                obj.OptimizationHistory(pos).Responses(iResponse).Value = r(iResponse);
            end
            obj.OptimizationHistory(pos).Objectives = obj.OptimizationObjectives;
            for iObjective = 1:obj.NumberOptimizationObjectives
                ind_obj = vertcat(obj.OptimizationResponses.Name) == obj.OptimizationObjectives(iObjective).Name;
                obj.OptimizationHistory(pos).Objectives(iObjective).Value = r(ind_obj);
            end
            obj.OptimizationHistory(pos).Constraints = obj.OptimizationConstraints;
            for iConstraint = 1:obj.NumberOptimizationConstraints
                ind_constr = vertcat(obj.OptimizationResponses.Name) == obj.OptimizationConstraints(iConstraint).Name;
                obj.OptimizationHistory(pos).Constraints(iConstraint).Value = r(ind_constr);
            end
            obj.OptimizationHistory(pos).Parameters = obj.OptimizationParameters;
            for iParameter = 1:obj.NumberOptimizationParameters
                obj.OptimizationHistory(pos).Parameters(iParameter).Value = x(iParameter);
            end
            obj.OptimizationHistory(pos).TimeElapsed = toc(obj.StartTime);
            OptimizationHistory = obj.OptimizationHistory;
            save(obj.SaveDirectory + "/OptimizationHistory.mat", 'OptimizationHistory'); 
        end

        function stop = plotFunction(obj, xnorm, optimValues, state)
            % TBD does scaling make sense here? check if iEl is OPT_System and not scale there?
            x = xnorm.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    obj.Elements(iEl).Element.plotFunction(x_iel, optimValues, state)
                end
            end
            stop = false;
        end

        % Resetting
        function resetOptimizationObjectives(obj)
            obj.OptimizationObjectives = struct('Name', {}, 'Index', {}, 'Scaling', {});
            obj.NumberOptimizationObjectives = 0;
        end
        
        function resetOptimizationConstraints(obj)
            obj.OptimizationConstraints = struct('Name', {}, 'Index', {}, 'Type', {}, 'Value', {}, 'Scaling', {});
            obj.NumberOptimizationConstraints = 0;
        end

        function resetOptimizationParameters(obj)
            resetOptimizationParameters@OPT_Element(obj);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    obj.Elements(iEl).Element.resetOptimizationParameters();
                end
            end
        end

        function resetOptimizationResponses(obj)
            resetOptimizationResponses@OPT_Element(obj);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    obj.Elements(iEl).Element.resetOptimizationResponses();
                end
            end
        end

        % Tests for responses
        % Numerical tests for derivatives: responses
        function [errAbs, errRel, r0, gradrAna, gradrNum] = checkResponses(obj, x0, numericalStep, testIndices)
            if ~exist("x0", "var") || isempty(x0)
                x0 = vertcat(obj.OptimizationParameters.Init);
            end
            if ~exist("numericalStep", "var") || isempty(numericalStep)
                numericalStep = (obj.UpperBounds - obj.LowerBounds)*1e-4;
            end
            if ~exist("testIndices", "var") || isempty(testIndices)
                testIndices = 1:obj.NumberOptimizationParameters;
            end
            % analytic responses and gradients
            [r0, gradrAna] = obj.response_function(x0);
            gradrAna = full(gradrAna);
            gradrAna = gradrAna(:, testIndices);
            % numerical gradients
            gradrNum = zeros(size(gradrAna));
            for iParam = 1:numel(testIndices)
                Step = numericalStep(iParam);
                param = testIndices(iParam);
                x1 = x0;
                x1(param) = x1(param) + Step;
            
                r1 = obj.response_function(x1);
                gradrNum(:, iParam) = (r1-r0)/Step;
                iParam
            end
            errAbs = abs(gradrNum-gradrAna);
            errRel = abs(gradrNum-gradrAna)./abs(gradrAna);
        end
    end 
end