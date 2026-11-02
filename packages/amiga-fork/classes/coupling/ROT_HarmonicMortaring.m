classdef ROT_HarmonicMortaring < CPL_HarmonicMortaring
    properties (Access = public)
        LoadFunction;
        Inertia;
        thetaDOF;
        omegaDOF;
        Beta = 1e-2;        % Friction coefficient
    end

    methods (Access = public)
        function obj = ROT_HarmonicMortaring(region1, sides1, region2, sides2, loadFunction, inertia)
            obj@CPL_HarmonicMortaring(region1, sides1, region2, sides2);
            obj.LoadFunction = loadFunction;
            obj.Inertia = inertia;
        end

        function setHarmonics(obj, type, args)
            arguments
                obj 
                type {mustBeMember(type, ["automatic", "manual"])};
                args.SinValues = [];
                args.CosValues = [];
                args.SymmetryType {mustBeMember(args.SymmetryType, ["periodic", "antiperiodic", "none"])} = "none";
                args.Symmetry = [];
                args.IncludeOffset (1,1) logical = false;
                args.DofRatio = 0.6;
            end
            cellargs = namedargs2cell(args);
            setHarmonics@CPL_HarmonicMortaring(obj, type, cellargs{:});
            % Angle and Velocity as additional DOFs
            obj.thetaDOF = obj.NumberDOFext+obj.NumberDOF + 1;
            obj.omegaDOF = obj.NumberDOFext+obj.NumberDOF + 2;
            obj.NumberDOF = obj.NumberDOF + 2;
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Mlin(obj.thetaDOF, obj.thetaDOF) = 1;
                Mlin(obj.omegaDOF, obj.omegaDOF) = obj.Inertia;
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            % Nonlinear because matrix can change with rotation angle
            if y(obj.thetaDOF) ~= obj.CouplingInterface2.RotationAngle
                obj.setRotationAngle(y(obj.thetaDOF));
            end
            Knonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            Knonlin(obj.a2DOFs, obj.couplingDOFs) = obj.Sgn2*obj.CouplingInterface2.CouplingMatrix *obj.MaterialParameter*obj.Region2.Length;
            Knonlin(obj.couplingDOFs, obj.a2DOFs) = obj.Sgn2*obj.CouplingInterface2.CouplingMatrix'*obj.MaterialParameter*obj.Region2.Length;
            if nargout == 2
                Klin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Klin(obj.a1DOFs, obj.couplingDOFs) = obj.Sgn1*obj.CouplingInterface1.CouplingMatrix *obj.MaterialParameter*obj.Region1.Length;
                Klin(obj.couplingDOFs, obj.a1DOFs) = obj.Sgn1*obj.CouplingInterface1.CouplingMatrix'*obj.MaterialParameter*obj.Region1.Length;
                Klin(obj.thetaDOF, obj.omegaDOF) = 1;
                Klin(obj.omegaDOF, obj.omegaDOF) = -obj.Beta;
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            if y(obj.thetaDOF) ~= obj.CouplingInterface2.RotationAngle
                obj.setRotationAngle(y(obj.thetaDOF));
            end
            Tem = -y(obj.a2DOFs)'*obj.CouplingInterface2.CouplingMatrixDer*y(obj.couplingDOFs)*obj.MaterialParameter*obj.Region2.Length;
            Fnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            Fnonlin(end) = Tem - obj.LoadFunction(t);
            if nargout == 2
                Flin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            end
        end

        function J = JacobianMatrix(obj, t, y)
            if y(obj.thetaDOF) ~= obj.CouplingInterface2.RotationAngle
                obj.setRotationAngle(y(obj.thetaDOF));
            end
            a2 = y(obj.a2DOFs);
            lambda = y(obj.couplingDOFs);
            Bstder = obj.Sgn2*obj.CouplingInterface2.CouplingMatrixDer*obj.MaterialParameter*obj.Region2.Length;
            Bstder2 = obj.Sgn2*obj.CouplingInterface2.CouplingMatrixDer2*obj.MaterialParameter*obj.Region2.Length;

            J = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            J(obj.a2DOFs, obj.thetaDOF) =       Bstder*lambda;
            J(obj.couplingDOFs, obj.thetaDOF) = Bstder'*a2;
            J(obj.omegaDOF, obj.a2DOFs) =       (Bstder*lambda)';
            J(obj.omegaDOF, obj.couplingDOFs) = a2'*Bstder;
            J(obj.omegaDOF, obj.thetaDOF) =     a2'*Bstder2*lambda;
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            a1 = ys(obj.a1DOFs, :);
            a2 = ys(obj.a2DOFs, :);
            obj.Solution.Values = ys(obj.couplingDOFs, :);
            obj.Solution.Angle = ys(obj.thetaDOF, :);
            obj.Region1.Solution.RotationAngle = obj.Solution.Angle;
            obj.Solution.Omega = ys(obj.omegaDOF, :);
            obj.Solution.Torque = zeros(size(ts));

            for iVal = 1:size(ys, 2)
                a1i = a1(:, iVal);
                a2i = a2(:, iVal);
                lambdai = obj.Solution.Values(:, iVal);
                obj.setRotationAngle(obj.Solution.Angle(iVal));
                obj.Solution.Torque(iVal) = -obj.MaterialParameter*a2i'*obj.CouplingInterface2.CouplingMatrixDer*lambdai*obj.Region2.Length;
            end
        end

        function plotSolution(obj)
            figure
            subplot(2,1,1);
            
            hold on
            plot(obj.Solution.Time, obj.Solution.Torque, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            xlabel("Time (s)");
            ylabel("Torque (Nm)");
            grid on
            
            subplot(2,1,2);
            hold on
            plot(obj.Solution.Time, obj.Solution.Omega, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            xlabel("Time (s)");
            ylabel("Angular Frequency (rad/s)");
            grid on
        end
    end
end