classdef CPL_HarmonicMortaring < EQS_Element
    properties
        MaterialParameter;
        Region1, Region2;
        CouplingInterface1, CouplingInterface2;
        a1DOFs, a2DOFs, couplingDOFs;
        RotationFunction;
        Sgn1, Sgn2;         % Sign that determins which region is outer/inner
    end

    methods
        function obj = CPL_HarmonicMortaring(region1, sides1, region2, sides2, args)
            arguments
                region1 {mustBeA(region1, {'IGA_Region', 'FEM_Region'})};
                sides1 (1,:) double;
                region2 {mustBeA(region2, {'IGA_Region', 'FEM_Region'})};
                sides2 (1,:) double;
                args.MaterialParameter = 1/(4*pi*1e-7);
                args.nGaussPoints = [];
                args.tolerance = 1e-10;
            end

            obj.MaterialParameter = args.MaterialParameter;
            if isa(region1, "IGA_Region")
                obj.CouplingInterface1 = IGA_CouplingInterface(region1, region1.getBoundaryPatchNumbers(sides1), args.nGaussPoints);
            else
                obj.CouplingInterface1 = FEM_CouplingInterface(region1, sides1);

            end
            if isa(region2, "IGA_Region")
                obj.CouplingInterface2 = IGA_CouplingInterface(region2, region2.getBoundaryPatchNumbers(sides2), args.nGaussPoints);
            else
                obj.CouplingInterface2 = FEM_CouplingInterface(region2, sides2);
            end
            
            if abs(obj.CouplingInterface1.Radius - obj.CouplingInterface2.Radius) > args.tolerance
                warning('CPL_HarmonicMortaring: Radii of coupling circles are different! \t Radius 1: %.10g \t Radius 2: %.10g \t Difference: %.10g \n', obj.CouplingInterface1.Radius, obj.CouplingInterface2.Radius, abs(obj.CouplingInterface1.Radius - obj.CouplingInterface2.Radius));
            end
            % check which region is outside and inside by comparing the first region
            if vecnorm(region1.getRegionCenter(1, true)) < vecnorm(region2.getRegionCenter(1, true))
                obj.Sgn1 = 1;
                obj.Sgn2 = -1;
            else
                obj.Sgn1 = -1;
                obj.Sgn2 = 1;
            end
            
            obj.Region1 = region1;
            obj.Region2 = region2;
            obj.addExternalElement(region1)
            if region1 ~= region2
                obj.addExternalElement(region2)
            else
                warning('CPL_HarmonicMortaring: You are coupling the same element to itself! This will results in wrong plots for rotated geometries!');
            end
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
            if type == "manual"
                sinValues = args.SinValues;
                cosValues = args.CosValues;
            elseif type == "automatic"
                if isempty(args.Symmetry)
                    args.Symmetry = round(2*pi*obj.CouplingInterface1.Radius/obj.CouplingInterface1.ArcLength);
                end
                maxDofs = args.DofRatio*min(obj.CouplingInterface1.NumberBoundaryDOFs, obj.CouplingInterface2.NumberBoundaryDOFs);
                nharm = maxDofs/2;
                if args.SymmetryType == "periodic"
                    harmonics = args.Symmetry:args.Symmetry:args.Symmetry*nharm;
                elseif args.SymmetryType == "antiperiodic"
                    harmonics = args.Symmetry/2:args.Symmetry:args.Symmetry*nharm;
                else
                    harmonics = 1:nharm;
                end
                sinValues = harmonics;
                if args.IncludeOffset
                    if args.SymmetryType == "antiperiodic"
                        warning("CPL_HarmonicMortaring: You specified an antiperiodic condition AND to include a possible offset at the interface. For antiperiodc conditions, this offset will always be zeros.")
                    end
                    cosValues = [0, harmonics];
                else
                    if args.SymmetryType == "periodic"
                        warning("CPL_HarmonicMortaring: You specified a periodic condition AND not to include a possible offset at the interface. The values at the interface will be zero in average, which might not be intended for periodic problems.")
                    end
                    cosValues = harmonics;
                end
            end
            obj.CouplingInterface1.setHarmonics(sinValues, cosValues);
            obj.CouplingInterface2.setHarmonics(sinValues, cosValues);
            obj.NumberDOF = numel(obj.CouplingInterface1.HarmonicsSin) + numel(obj.CouplingInterface1.HarmonicsCos);
            obj.a1DOFs = 1:obj.Region1.NumberDOF;
            obj.a2DOFs = obj.Region1.NumberDOF+(1:obj.Region2.NumberDOF);
            obj.couplingDOFs = obj.NumberDOFext+(1:obj.NumberDOF);
        end

        function setRotationAngle(obj, angleRad)
            obj.CouplingInterface2.setRotationAngle(angleRad);
            % rotate region
            obj.Region1.RotationAngle = angleRad;
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            if ~isempty(obj.RotationFunction)
                obj.setRotationAngle(obj.RotationFunction(t))
            end
            % Nonlinear because matrix can change with rotation angle
            Knonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            Knonlin(obj.a2DOFs, obj.couplingDOFs) = obj.Sgn2*obj.CouplingInterface2.CouplingMatrix * obj.MaterialParameter * obj.Region2.Length;
            Knonlin(obj.couplingDOFs, obj.a2DOFs) = obj.Sgn2*obj.CouplingInterface2.CouplingMatrix' * obj.MaterialParameter * obj.Region2.Length;
            if nargout == 2
                Klin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Klin(obj.a1DOFs, obj.couplingDOFs) = obj.Sgn1*obj.CouplingInterface1.CouplingMatrix * obj.MaterialParameter * obj.Region1.Length;
                Klin(obj.couplingDOFs, obj.a1DOFs) = obj.Sgn1*obj.CouplingInterface1.CouplingMatrix' * obj.MaterialParameter * obj.Region1.Length;
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            if nargout == 2
                Flin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            end
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
        end

        function postprocess(obj, ts, ys, ~)
            a1 = ys(obj.a1DOFs, :);
            a2 = ys(obj.a2DOFs, :);
            obj.Solution.Time = ts;
            obj.Solution.Values = ys(obj.couplingDOFs, :);
            obj.Solution.Torque = zeros(size(ts));
            obj.Region1.Solution.RotationAngle = zeros(1, numel(ts));
            for iT = 1:numel(ts)
                t = ts(iT);
                if ~isempty(obj.RotationFunction)
                    obj.setRotationAngle(obj.RotationFunction(t));
                    obj.Region1.Solution.RotationAngle(iT) = obj.RotationFunction(t);
                end
                obj.Solution.Torque(iT) = -obj.MaterialParameter*a2(:,iT)'*obj.CouplingInterface2.CouplingMatrixDer*obj.Solution.Values(:,iT)*obj.Region2.Length;
            end
        end

        function plotSolution(obj)
            figure()
            plot(obj.Solution.Time, obj.Solution.Torque, "Color", TUDa_getColor("1b"), LineWidth=1.5);
            ylabel("Torque (Nm)")
            xlabel("Time (s)")
            grid on
        end

        function plotCouplingValues(obj)
            figure;
            idxsin = mod(obj.CouplingInterface1.CouplingIndices, 2) == 1;
            idxcos = mod(obj.CouplingInterface1.CouplingIndices, 2) == 0;
            sinValues = obj.Solution.Values(idxsin);
            cosValues = obj.Solution.Values(idxcos);
            xlims = [min(obj.CouplingInterface1.HarmonicsAll)-1, max(obj.CouplingInterface1.HarmonicsAll)+1];
            subplot(2,1,1)
            bar(obj.CouplingInterface1.HarmonicsSin, sinValues);
            xticks(obj.CouplingInterface1.HarmonicsSin);
            xlim(xlims);
            title("Sin Values");
            subplot(2,1,2)
            bar(obj.CouplingInterface1.HarmonicsCos, cosValues);
            xticks(obj.CouplingInterface1.HarmonicsCos);
            xlim(xlims);
            title("Cos Values");
        end
    end
end