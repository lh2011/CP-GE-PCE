classdef CPL_HarmonicMortaring_OPT < OPT_Element & CPL_HarmonicMortaring

    properties
        RadiusName
    end

    methods
        function obj = CPL_HarmonicMortaring_OPT(region1, sides1, region2, sides2, args)
            arguments
                region1 {mustBeA(region1, {'IGA_Region', 'FEM_Region'})};
                sides1 (1,:) double;
                region2 {mustBeA(region2, {'IGA_Region', 'FEM_Region'})};
                sides2 (1,:) double;
                args.MaterialParameter = 1/(4*pi*1e-7);
                args.nGaussPoints = [];
                args.tolerance = 1e-10;
            end
            unpackedargs = namedargs2cell(args);
            obj@CPL_HarmonicMortaring(region1, sides1, region2, sides2, unpackedargs{:});
        end

        function optimizationUpdate(obj, x)
            obj.CouplingInterface1.generateSpacesAndMeshes();
            obj.CouplingInterface2.generateSpacesAndMeshes();
            obj.CouplingInterface1.generateCouplingMatrices();
            obj.CouplingInterface2.generateCouplingMatrices();

            indAngle = findStructName(obj.OptimizationParameters, "RotationAngle");
            if ~isempty(indAngle)
                obj.setRotationAngle(x(indAngle));
            end
        end

        function r = responses(obj, t, x, y)
            r = zeros(obj.NumberOptimizationResponses, 1);

            a2 = y(obj.a2DOFs);
            lambda = y(obj.couplingDOFs);

            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            if ~isempty(indTorque)
                nu0 = obj.MaterialParameter;
                r(indTorque) = obj.Sgn2*nu0*a2'*obj.CouplingInterface2.CouplingMatrixDer*lambda*obj.Region2.Length;
            end
        end

        function drdy = partial_r_partial_y(obj, t, x, y)
            drdy = sparse(obj.NumberOptimizationResponses, obj.NumberDOF + obj.NumberDOFext);

            a2 = y(obj.a2DOFs);
            lambda = y(obj.couplingDOFs);

            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            if ~isempty(indTorque)
                nu0 = obj.MaterialParameter;
                drdy(indTorque, :) = obj.Sgn2*nu0*[sparse(obj.Region1.NumberDOF, 1); ...
                    obj.CouplingInterface2.CouplingMatrixDer*lambda; ...
                    obj.CouplingInterface2.CouplingMatrixDer'*a2] *obj.Region2.Length;
            end
        end

        function Lx = function_Lx(obj, t, x, y, z)
            Lx = zeros(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);
            indAngle = findStructName(obj.OptimizationParameters, "RotationAngle");
            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            if  ~isempty(indAngle) && ~isempty(indTorque)
                z_Torque = z(:,indTorque);
                a2 = y(obj.a2DOFs);
                lambda = y(obj.couplingDOFs);
                Lx(indTorque, indAngle) = z_Torque(obj.a2DOFs)'*obj.Sgn2*obj.CouplingInterface2.CouplingMatrixDer * obj.MaterialParameter * obj.Region2.Length * lambda + ...
                                          z_Torque(obj.couplingDOFs)'*obj.Sgn2*obj.CouplingInterface2.CouplingMatrixDer' * obj.MaterialParameter * obj.Region2.Length * a2;
            end

        end

        function drdx = partial_r_partial_x(obj, t, x, y)
            drdx = sparse(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);

            indTorque = findStructName(obj.OptimizationResponses, "Torque");
            indAngle = findStructName(obj.OptimizationParameters, "RotationAngle");
            if ~isempty(indTorque) && ~isempty(indAngle)
                a2 = y(obj.a2DOFs);
                lambda = y(obj.couplingDOFs);
                nu0 = obj.MaterialParameter;
                drdx(indTorque,indAngle) = obj.Sgn2*nu0*a2'*obj.CouplingInterface2.CouplingMatrixDer2*lambda*obj.Region2.Length;
            end
            indRadius = findStructName(obj.OptimizationParameters, obj.RadiusName);
            if ~isempty(indTorque) && ~isempty(indRadius)
                a2 = y(obj.a2DOFs);
                lambda = y(obj.couplingDOFs);
                nu0 = obj.MaterialParameter;
                R = obj.CouplingInterface2.Radius;
                drdx(indTorque,indRadius) = 1/R*obj.Sgn2*nu0*a2'*obj.CouplingInterface2.CouplingMatrixDer*lambda*obj.Region2.Length;
            end
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            dedx = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberOptimizationParameters);
            indAngle = findStructName(obj.OptimizationParameters, "RotationAngle");
            if  ~isempty(indAngle)
                a2 = y(obj.a2DOFs);
                lambda = y(obj.couplingDOFs);
                dedx(obj.a2DOFs, indAngle) = obj.Sgn2*obj.CouplingInterface2.CouplingMatrixDer * obj.MaterialParameter * obj.Region2.Length * lambda;
                dedx(obj.couplingDOFs, indAngle) = obj.Sgn2*obj.CouplingInterface2.CouplingMatrixDer' * obj.MaterialParameter * obj.Region2.Length * a2;
            end

            indRadius = findStructName(obj.OptimizationParameters, obj.RadiusName);
            if  ~isempty(indRadius)
                a1 = y(obj.a1DOFs);
                a2 = y(obj.a2DOFs);
                lambda = y(obj.couplingDOFs);
                R = obj.CouplingInterface1.Radius;
                dedx(obj.a2DOFs, indRadius) = 1/R*obj.Sgn2*obj.CouplingInterface2.CouplingMatrix * obj.MaterialParameter * obj.Region2.Length * lambda;
                dedx(obj.couplingDOFs, indRadius) = 1/R*obj.Sgn2*obj.CouplingInterface2.CouplingMatrix' * obj.MaterialParameter * obj.Region2.Length * a2;
                dedx(obj.a1DOFs, indRadius) = 1/R*obj.Sgn1*obj.CouplingInterface1.CouplingMatrix * obj.MaterialParameter * obj.Region1.Length * lambda;
                dedx(obj.couplingDOFs, indRadius) = dedx(obj.couplingDOFs, indRadius) + 1/R*obj.Sgn1*obj.CouplingInterface1.CouplingMatrix' * obj.MaterialParameter * obj.Region1.Length * a1;
            end
        end

        function addOptimizationRadius(obj, name, initVal, minVal, maxVal)
            obj.addOptimizationParameter(name, initVal, minVal, maxVal);
            obj.RadiusName = name;
        end

        function plotFunction(obj, x, optimValues, state)
            % do nothing
        end
    end
end

