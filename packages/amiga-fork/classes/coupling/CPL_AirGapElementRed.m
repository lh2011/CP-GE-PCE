% reduced variant for CPL_AirGapElement where the internal DOFs are
% eliminated, which makes it symmetric.

classdef CPL_AirGapElementRed < CPL_AirGapElement
    properties
        KFinv;
        NumberInt;  % replaces NumberDOF, as they are removed
    end

    methods
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
            setHarmonics@CPL_AirGapElement(obj, type, cellargs{:});
            obj.NumberInt = obj.NumberDOF;
            obj.NumberDOF = 0;
            obj.KFinv = inv(obj.KF);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            if ~isempty(obj.RotationFunction)
                obj.setRotationAngle(obj.RotationFunction(t))
            end
            Klin = sparse(obj.NumberDOFext, obj.NumberDOFext);
            Knonlin = sparse(obj.NumberDOFext, obj.NumberDOFext);

            B1 = obj.CouplingInterface1.CouplingMatrix*obj.Region1.Length;
            B2 = obj.CouplingInterface2.CouplingMatrix*obj.Region2.Length;
            Knonlin(obj.a1DOFs, obj.a1DOFs) = B1*obj.C1*(obj.KF\(obj.D1*B1'));
            Knonlin(obj.a1DOFs, obj.a2DOFs) = B1*obj.C1*(obj.KF\(obj.D2*B2'));
            Knonlin(obj.a2DOFs, obj.a1DOFs) = B2*obj.C2*(obj.KF\(obj.D1*B1'));
            Knonlin(obj.a2DOFs, obj.a2DOFs) = B2*obj.C2*(obj.KF\(obj.D2*B2'));
        end

        function postprocess(obj, ts, ys, ~)
            obj.Solution.Time = ts;
            obj.Solution.Values = zeros(obj.NumberInt, numel(ts));
            art = ys(obj.a1DOFs, :);
            ast = ys(obj.a2DOFs, :);
            for iValue = 1:numel(ts)
                t = ts(iValue);
                if ~isempty(obj.RotationFunction)
                    obj.setRotationAngle(obj.RotationFunction(t))
                    obj.Region1.Solution.RotationAngle(iValue) = obj.RotationFunction(t);
                end
                B1 = obj.CouplingInterface1.CouplingMatrix*obj.Region1.Length;
                B2 = obj.CouplingInterface2.CouplingMatrix*obj.Region2.Length;
                obj.Solution.Values(:, iValue) = obj.KF\(obj.D1*B1'*art(:,iValue) + obj.D2*B2'*ast(:,iValue));
            end
            obj.Solution.Torque = obj.calcTorqueCoupling(obj.Solution.Values);
        end
    end
end