classdef  OPT_System_2nd < OPT_System & OPT_Element_2nd
    methods
        % Overwrite optimize with Hessian function for options
        function x_opt = optimize(obj, x0, newOptions)
            if ~exist("x0", "var") || isempty(x0)
                x0 = vertcat(obj.OptimizationParameters.Init);
            end
            % standard options for second order optimization
            standardOptions.HessianFcn = @(x, lambda) obj.hessian_objectives(x, lambda);
            % standardOptions.InitBarrierParam = 1e-6;
            % standarsOptions.BarrierParamUpdate = 'predictor-corrector';

            if exist("newOptions", "var")
                field_names = fieldnames(newOptions);
                for iNewOpts = 1:numel (field_names)
                    standardOptions.(field_names{iNewOpts}) = newOptions.(field_names{iNewOpts});
                end
            end

            % '', 1e-6, ...
            % 'StepTolerance', 1e-8, ... % 
            %     'OptimalityTolerance', 1e-6, ... 
            %     'ConstraintTolerance', 1e-6, ... 

            x_opt = optimize@OPT_System(obj, x0, standardOptions);
        end

        function H = hessian_objectives(obj, x_norm, lambda, t)
            if ~exist("t", "var")
                t = 0;
            end
            % Rescale values from 0-1 to real values
            x = x_norm.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            H = zeros(obj.NumberOptimizationParameters, obj.NumberOptimizationParameters);

            Hr = obj.hessian_responses(x, t);
            for iObjective = 1:numel(obj.OptimizationObjectives)
                ind_obj = obj.OptimizationObjectives(iObjective).Index; 
                scale = obj.OptimizationObjectives(iObjective).Scaling; 
                H = H + squeeze(Hr(ind_obj, :, :))*scale;
            end
            
            iIneq = 0;
            iEq = 0;
            for iConstraint = 1:numel(obj.OptimizationConstraints)
                ind_constr = obj.OptimizationConstraints(iConstraint).Index; 
                scale = obj.OptimizationConstraints(iConstraint).Scaling;
                type = obj.OptimizationConstraints(iConstraint).Type;

                if type == "="
                    iEq = iEq + 1;
                    H = H + lambda.eqnonlin(iEq)*squeeze(Hr(ind_constr,:,:))*scale;
                elseif type == "<"
                    iIneq = iIneq + 1;
                    H = H + lambda.ineqnonlin(iIneq)*squeeze(Hr(ind_constr,:,:))*scale;
                elseif type == ">"
                    iIneq = iIneq + 1;
                    H = H - lambda.ineqnonlin(iIneq)*squeeze(Hr(ind_constr,:,:))*scale;
                else
                    error("wrong type");
                end
            end

            H = H .* ((obj.UpperBounds-obj.LowerBounds)*(obj.UpperBounds-obj.LowerBounds)');
        end

        function [Hr, gradr, r] = hessian_responses(obj, x, t)
            if ~exist("t", "var")
                t = 0;
            end
            mytimer("SecondOrder", "start")
            mytimer("SecondOrder-Build", "start")
            obj.optimizationUpdate(x);
            mytimer("SecondOrder-Build", "stop")
            
            if isempty(obj.lasty)
                y0 = zeros(obj.NumberDOF, 1);
            else
                y0 = obj.lasty;
            end
            mytimer("SecondOrder-Newton", "start")
            [y, Jac] = obj.solveStatic(t, y0, "Verbose", false);
            mytimer("SecondOrder-Newton", "stop")
            obj.lasty = y;

            mytimer("SecondOrder-Gradient", "start")
            dedx = obj.partial_e_partial_x(t, x, y);
            drdy = obj.partial_r_partial_y(t, x, y);
            % Applying boundary conditions
            dedxBC = dedx(obj.K, :) + obj.G'*(dedx(obj.R, :));
            drdyBC = drdy(:, obj.K) + drdy(:, obj.R)*obj.G;
            
            % Jac_dec = decomposition(Jac, "LDL", "upper");
            Jac_dec = decomposition(Jac, "LU");

            % zk = -(Jac'\drdyBC');
            zk = -(Jac_dec'\drdyBC');
            zr = obj.G*zk;
            z = zeros(obj.NumberDOF, obj.NumberOptimizationResponses);
            z(obj.K, :) = zk;
            z(obj.R, :) = zr;

            % If first derivative is needed:
            if nargout > 1
                r = obj.responses(t, x, y);
                drdx = obj.partial_r_partial_x(t, x, y);
                gradr = drdx + zk'*dedxBC; 
            end
            mytimer("SecondOrder-Gradient", "stop")

            % Second derivatives
            mytimer("SecondOrder-Lxx", "start")
            Lxx = obj.function_Lxx(t, x, y, z);
            mytimer("SecondOrder-Lxx", "stop")

            mytimer("SecondOrder-Lyy", "start")
            Lyy = obj.function_Lyy(t, x, y, z);
            LyyBC = Lyy(obj.K, obj.K, :) + ...                
                permute(ttt(Lyy(obj.K, obj.R, :), sptensor(obj.G), [2], [1]), [1,3,2]) + ...                
                permute(ttt(Lyy(obj.R, obj.K, :), sptensor(obj.G), [1], [1]), [3,1,2]) + ...
                permute(ttt(ttt(Lyy(obj.R, obj.R,:), sptensor(obj.G), 1, 1), sptensor(obj.G), 1, 1), [2,3,1]);
            mytimer("SecondOrder-Lyy", "stop")

            mytimer("SecondOrder-Lxy", "start")
            Lxy = obj.function_Lxy(t, x, y, z);

            LxyBC = Lxy(:, obj.K, :) + permute(ttt(Lxy(:, obj.R, :), sptensor(obj.G) , [2], [1]), [1,3,2]);

            mytimer("SecondOrder-Lxy", "stop")

            mytimer("SecondOrder-Hessian", "start")
            % Tangent equation
            dydx = Jac_dec\(-dedxBC);       % direct method
            Hr = zeros(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters, obj.NumberOptimizationParameters);
            for iResponse = 1:numel(obj.OptimizationResponses)
                Lyyi = double(sptenmat(LyyBC(:, :, iResponse), 1));
                Lxyi = double(sptenmat(LxyBC(:, :, iResponse), 1));
                Lxxi = Lxx(:, :, iResponse);
                % Auxiliary adjoint equation
                dzdx = -Jac_dec'\(Lxyi' + Lyyi*dydx);
                Hr(iResponse, :, :) = Lxxi + Lxyi*dydx + dedxBC'*dzdx;
            end
            mytimer("SecondOrder-Hessian", "stop")
            mytimer("SecondOrder", "stop")
        end
      
        function Lyy = function_Lyy(obj, t, x, y, z)
            Lyy = sptensor([], [], [obj.NumberDOF, obj.NumberDOF, obj.NumberOptimizationResponses]);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element_2nd")
                    ind_y = obj.Elements(iEl).Indices;
                    y_iel = y(ind_y);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    ind_r = obj.Elements(iEl).IndicesResponses;
                    z_iel = z(ind_y, ind_r);
                    
                    Lyy_iel = obj.Elements(iEl).Element.function_Lyy(t, x_iel, y_iel, z_iel);
                    if ~isempty(Lyy_iel) && nnz(Lyy_iel) > 0
                        Lyy(ind_y, ind_y, ind_r) = Lyy(ind_y, ind_y, ind_r) + squeeze(Lyy_iel);
                    end
                end
            end
        end

        function Lxx = function_Lxx(obj, t, x, y, z)
            Lxx = zeros([obj.NumberOptimizationParameters, obj.NumberOptimizationParameters, obj.NumberOptimizationResponses]);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element_2nd")
                    ind_y = obj.Elements(iEl).Indices;
                    y_iel = y(ind_y);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    ind_r = obj.Elements(iEl).IndicesResponses;
                    z_iel = z(ind_y, ind_r);
                    
                    Lxx_iel = obj.Elements(iEl).Element.function_Lxx(t, x_iel, y_iel, z_iel);
                    if ~isempty(Lxx_iel)
                        Lxx(ind_x, ind_x, ind_r) = Lxx(ind_x, ind_x, ind_r) + squeeze(Lxx_iel);
                    end
                end
            end
        end

        function Lxy = function_Lxy(obj, t, x, y, z)
            Lxy = sptensor([], [], [obj.NumberOptimizationParameters, obj.NumberDOF, obj.NumberOptimizationResponses]);
            % Lxy = zeros(obj.NumberOptimizationParameters, obj.NumberDOF, obj.NumberOptimizationResponses);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element_2nd")
                    ind_y = obj.Elements(iEl).Indices;
                    y_iel = y(ind_y);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    ind_r = obj.Elements(iEl).IndicesResponses;
                    z_iel = z(ind_y, ind_r);
                    
                    Lxy_iel = obj.Elements(iEl).Element.function_Lxy(t, x_iel, y_iel, z_iel);
                    if ~isempty(Lxy_iel)
                        Lxy(ind_x, ind_y, ind_r) = Lxy(ind_x, ind_y, ind_r) + squeeze(Lxy_iel);
                        % Lxy(ind_x, ind_y, ind_r) = Lxy(ind_x, ind_y, ind_r) + Lxy_iel;
                    end
                end
            end
        end
    end 
end