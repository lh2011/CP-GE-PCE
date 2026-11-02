classdef CPL_AirGapElement < EQS_Element
    properties
        MaterialParameter;
        Region1, Region2;
        CouplingInterface1, CouplingInterface2;
        HarmonicsSin, HarmonicsCos;
        a1DOFs;
        a2DOFs;
        internalDOFs;
        ri, ro;

        RotationFunction;
        % scale=1;
        AirGapThickness;
        C1, C2;
        D1, D2;
        KF;
    end

    methods
        function obj = CPL_AirGapElement(region1, sides1, region2, sides2, args)
            arguments
                region1 {mustBeA(region1, {'IGA_Region', 'FEM_Region'})};
                sides1 (1,:) double;
                region2 {mustBeA(region2, {'IGA_Region', 'FEM_Region'})};
                sides2 (1,:) double;
                args.MaterialParameter = 1/(4*pi*1e-7);
                args.nGaussPoints = [];
                args.tolerance = 1e-10;
            end
            assert(region1.Length == region2.Length, 'CPL_AirGapElement: The Air Gap Element requires regions 1 and 2 to be of equal length!');
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
            
            obj.AirGapThickness = abs(obj.CouplingInterface1.Radius - obj.CouplingInterface2.Radius);
            if obj.AirGapThickness < args.tolerance
                warning('CPL_AirGapElement: Radii of coupling circles are almost the same! \t Radius 1: %.10g \t Radius 2: %.10g \t Air Gap %.10g \n', obj.CouplingInterface1.Radius, obj.CouplingInterface2.Radius, obj.AirGapThickness);
            end

            % conditioning
            obj.ri = min(obj.CouplingInterface1.Radius, obj.CouplingInterface2.Radius);
            obj.ro = max(obj.CouplingInterface1.Radius, obj.CouplingInterface2.Radius);
                
            obj.Region1 = region1;
            obj.Region2 = region2;
            obj.addExternalElement(region1)
            if region1 ~= region2
                obj.addExternalElement(region2)
            else
                warning('CPL_AirGapElement: You are coupling the same element to itself! This will results in wrong plots for rotated geometries!');
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
                        warning("CPL_AirGapElement: You specified an antiperiodic condition AND to include a possible offset at the interface. For antiperiodc conditions, this offset will always be zeros.")
                    end
                    cosValues = [0, harmonics];
                else
                    if args.SymmetryType == "periodic"
                        warning("CPL_AirGapElement: You specified a periodic condition AND not to include a possible offset at the interface. The values at the interface will be zero in average, which might not be intended for periodic problems.")
                    end
                    cosValues = harmonics;
                end
            end
            obj.CouplingInterface1.setHarmonics(sinValues, cosValues);
            obj.CouplingInterface2.setHarmonics(sinValues, cosValues);
            obj.HarmonicsSin = obj.CouplingInterface1.HarmonicsSin;
            obj.HarmonicsCos = obj.CouplingInterface1.HarmonicsCos;
            obj.NumberDOF = 2*(numel(obj.HarmonicsSin) + numel(obj.HarmonicsCos));
            obj.a1DOFs = 1:obj.Region1.NumberDOF;
            obj.a2DOFs = obj.Region1.NumberDOF+(1:obj.Region2.NumberDOF);
            obj.internalDOFs = obj.NumberDOFext+(1:obj.NumberDOF);
            obj.createMatrices();
        end

        function createMatrices(obj)
            ns = obj.CouplingInterface1.HarmonicsAll;
            n = numel(ns);
            usedRows = obj.CouplingInterface1.CouplingIndices;
            usedCols = union(2*usedRows-1, 2*usedRows);

            % Factor matrix C1 for derivatives at interface 1
            R1 = obj.CouplingInterface1.Radius;
            sin_fact1 = ns./R1.*(R1/obj.ro).^ns;
            sin_fact2 = -ns./R1.*(obj.ri/R1).^ns;
            cos_fact1 = ns./R1.*(R1/obj.ro).^ns;
            cos_fact2 = -ns./R1.*(obj.ri/R1).^ns;
            % zero mode correction
            if ns(1) == 0
                cos_fact1(1) = 0;
                cos_fact2(1) = 1/R1;
            end
            factors = reshape([sin_fact1; sin_fact2; cos_fact1; cos_fact2], 1, []);
            rows = reshape(repmat(1:2*n, 2,1), 1, []);
            cols = 1:4*n;

            if obj.CouplingInterface2.Radius > obj.CouplingInterface1.Radius
                sgn1 = 1;
                sgn2 = -1;
            else
                sgn1 = -1;
                sgn2 = 1;
            end

            obj.C1 = sparse(rows, cols, factors, 2*n, 4*n);
            obj.C1 = sgn1*obj.C1(usedRows, usedCols)*obj.MaterialParameter;

            % Factor matrix C2 for derivatives at interface 2
            R2 = obj.CouplingInterface2.Radius;
            sin_fact1 = ns./R2.*(R2/obj.ro).^ns;
            sin_fact2 = -ns./R2.*(obj.ri/R2).^ns;
            cos_fact1 = ns./R2.*(R2/obj.ro).^ns;
            cos_fact2 = -ns./R2.*(obj.ri/R2).^ns;
            % zero mode correction
            if ns(1) == 0
                cos_fact1(1) = 0;
                cos_fact2(1) = 1/R2;
            end
            factors = reshape([sin_fact1; sin_fact2; cos_fact1; cos_fact2], 1, []);
            rows = reshape(repmat(1:2*n, 2,1), 1, []);
            cols = 1:4*n;

            obj.C2 = sparse(rows, cols, factors, 2*n, 4*n);
            obj.C2 = sgn2*obj.C2(usedRows, usedCols)*obj.MaterialParameter;
            
            % Factor matrix KF for testing analytic function with sin/cos
            L1 = obj.CouplingInterface1.ArcLength;
            sin_fact1 = 0.5*L1*(R1/obj.ro).^ns;
            sin_fact2 = 0.5*L1*(obj.ri/R1).^ns;
            cos_fact1 = 0.5*L1*(R1/obj.ro).^ns;
            cos_fact2 = 0.5*L1*(obj.ri/R1).^ns;
            % zero mode correction
            if ns(1) == 0
                cos_fact1(1) = L1;
                cos_fact2(1) = L1*log(R1);
            end
            factors = reshape([sin_fact1; sin_fact2; cos_fact1; cos_fact2], 1, []);
            factorR1 = factors(usedCols) * obj.Region1.Length;

            L2 = obj.CouplingInterface2.ArcLength;
            sin_fact1 = 0.5*L2*(R2/obj.ro).^ns;
            sin_fact2 = 0.5*L2*(obj.ri/R2).^ns;
            cos_fact1 = 0.5*L2*(R2/obj.ro).^ns;
            cos_fact2 = 0.5*L2*(obj.ri/R2).^ns;
            % zero mode correction
            if ns(1) == 0
                cos_fact1(1) = L2;
                cos_fact2(1) = L2*log(R2);
            end
            factors = reshape([sin_fact1; sin_fact2; cos_fact1; cos_fact2], 1, []);
            factorR2 = factors(usedCols) * obj.Region2.Length;

            obj.KF = sparse(reshape(repmat(1:2:obj.NumberDOF, [2,1]), 1, []), [1:obj.NumberDOF], factorR1, obj.NumberDOF, obj.NumberDOF) + ...
                     sparse(reshape(repmat(2:2:obj.NumberDOF, [2,1]), 1, []), [1:obj.NumberDOF], factorR2, obj.NumberDOF, obj.NumberDOF);

            % Factor matrix for D1 u_num = u_ana
            obj.D1 = sparse(1:2:4*n, 1:2*n, ones(2*n,1), 4*n,2*n);
            obj.D1 = obj.D1(usedCols, usedRows);
            obj.D2 = sparse(2:2:4*n, 1:2*n, ones(2*n,1), 4*n,2*n);
            obj.D2 = obj.D2(usedCols, usedRows);
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
            
            % CouplingMatrix 2 may change with rotation angle
            Knonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            Knonlin(obj.a2DOFs, obj.internalDOFs) = obj.CouplingInterface2.CouplingMatrix*obj.C2*obj.Region2.Length;
            Knonlin(obj.internalDOFs, obj.a2DOFs) = obj.D2*obj.CouplingInterface2.CouplingMatrix'*obj.Region2.Length;

            if nargout == 2
                Klin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Klin(obj.a1DOFs, obj.internalDOFs) = obj.CouplingInterface1.CouplingMatrix*obj.C1*obj.Region1.Length;
                Klin(obj.internalDOFs, obj.a1DOFs) = obj.D1*obj.CouplingInterface1.CouplingMatrix'*obj.Region1.Length;
                Klin(obj.internalDOFs, obj.internalDOFs) = -obj.KF;
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
            obj.Solution.Time = ts;
            obj.Solution.Values = ys(obj.internalDOFs, :);
            obj.Solution.Torque = obj.calcTorqueCoupling(ys(obj.internalDOFs, :));
            obj.Region1.Solution.RotationAngle = zeros(1, numel(ts));
            for iValue = 1:numel(ts)
                t = ts(iValue);
                if ~isempty(obj.RotationFunction)
                    obj.Region1.Solution.RotationAngle(iValue) = obj.RotationFunction(t);
                end
            end
        end

        function T = calcTorqueBrBtRt(obj)
            T = 0;
            a = obj.Region1.Solution.Values;
            for iBnd = 1:numel(obj.CouplingInterface1.PatchBoundaries)
                patchDOFs = obj.Region1.Spaces.gnum{obj.Region1.PatchBoundaries(obj.CouplingInterface1.PatchBoundaries(iBnd)).patches};
                A_z = a(patchDOFs);
                T = T + obj.MaterialParameter * obj.Region1.Length * ...
                op_Br_Bt_dGamma(obj.CouplingInterface1.SpacesEval{iBnd}, obj.CouplingInterface1.MeshesEval{iBnd}, A_z);
            end
        end

        function T = calcTorqueBrBtSt(obj)
            T = 0;
            a = obj.Region2.Solution.Values;
            for iBnd = 1:numel(obj.CouplingInterface2.PatchBoundaries)
                patchDOFs = obj.Region2.Spaces.gnum{obj.Region2.PatchBoundaries(obj.CouplingInterface2.PatchBoundaries(iBnd)).patches};
                A_z = a(patchDOFs);
                T = T + obj.MaterialParameter * obj.Region2.Length * ...
                op_Br_Bt_dGamma(obj.CouplingInterface2.SpacesEval{iBnd}, obj.CouplingInterface2.MeshesEval{iBnd}, A_z);
            end
        end

        function Torques = calcTorqueCoupling(obj, ys)

            TorqueMatrix = sparse(obj.NumberDOF, obj.NumberDOF);
            commonSinCos = intersect(obj.HarmonicsSin, obj.HarmonicsCos);
            allSinCos = sort([obj.HarmonicsSin,obj.HarmonicsSin,obj.HarmonicsCos,obj.HarmonicsCos]);
            % assert(numel(allSinCos) == obj.NumberDOF);

            [~, indexSin, ~] = intersect(allSinCos, commonSinCos);

            % find harmonics where both sin and cos values are used
            for iHarm = 1:numel(commonSinCos)
                n = commonSinCos(iHarm);
                idx = indexSin(iHarm):indexSin(iHarm)+3;
                TorqueMatrix(idx, idx) = (obj.ri/obj.ro)^n* ...
                                          [0,  0,  0,  n^2;
                                           0,  0, -n^2,0;
                                           0, -n^2,0,  0;
                                           n^2,0,  0,  0];
            end
            Torques = zeros(1, size(ys, 2));
            theta = obj.CouplingInterface1.ArcLength/obj.CouplingInterface1.Radius;
            for iValues = 1:size(ys, 2)
                ys_i = ys(:, iValues);
                Lz = (obj.Region1.Length + obj.Region2.Length)/2;
                Torques(iValues) = 0.5*theta*Lz*obj.MaterialParameter*ys_i'*TorqueMatrix*ys_i;
            end            
        end

        function plotCouplingValues(obj)
            figure;

            CouplingIndices = obj.CouplingInterface1.CouplingIndices;
            used_indices = union(2*CouplingIndices, 2*CouplingIndices-1);

            sinValues1 = obj.Solution.Values((mod(used_indices, 4) == 1));
            sinValues2 = obj.Solution.Values((mod(used_indices, 4) == 2));
            cosValues1 = obj.Solution.Values((mod(used_indices, 4) == 3));
            cosValues2 = obj.Solution.Values((mod(used_indices, 4) == 0));
            xlims = [min(obj.CouplingInterface1.HarmonicsAll)-1, max(obj.CouplingInterface1.HarmonicsAll)+1];
            subplot(4,1,1)
            bar(obj.HarmonicsSin, sinValues1);
            xticks(obj.HarmonicsSin);
            xlim(xlims);
            title("Sin Values 1");
            subplot(4,1,2)
            bar(obj.HarmonicsSin, sinValues2);
            xticks(obj.HarmonicsSin);
            xlim(xlims);
            title("Sin Values 2");

            subplot(4,1,3)
            bar(obj.HarmonicsCos, cosValues1);
            xticks(obj.HarmonicsCos);
            xlim(xlims);
            title("Cos Values 1");
            subplot(4,1,4)
            bar(obj.HarmonicsCos, cosValues2);
            xticks(obj.HarmonicsCos);
            xlim(xlims);
            title("Cos Values 2");
        end

        function plotSolution(obj, interval)
            arguments
                obj 
                interval = [0, 2*pi];
            end
            ns = obj.CouplingInterface1.HarmonicsAll;
            n = numel(ns);
            allFunctions = cell(4, n);
            for iFunc = 1:n
                ni = ns(iFunc);
                allFunctions{1, iFunc} = @(r, theta) sin(ni*theta).*(r./obj.ro).^ni;
                allFunctions{2, iFunc} = @(r, theta) sin(ni*theta).*(obj.ri./r).^(ni);
                allFunctions{3, iFunc} = @(r, theta) cos(ni*theta).*(r./obj.ro).^(ni);
                allFunctions{4, iFunc} = @(r, theta) cos(ni*theta).*(obj.ri./r).^(ni);
            end
            if ns(1) == 0
                allFunctions{1,1} = @(r, theta) zeros(size(r));
                allFunctions{2,1} = @(r, theta) zeros(size(r));
                allFunctions{3,1} = @(r, theta) ones(size(r));
                allFunctions{4,1} = @(r, theta) log(r);
            end
            allFunctions = reshape(allFunctions, 1, []);
            CouplingIndices = obj.CouplingInterface1.CouplingIndices;
            used_indices = union(2*CouplingIndices, 2*CouplingIndices-1);

            r = linspace(obj.CouplingInterface1.Radius, obj.CouplingInterface2.Radius, 10);
            thetas = linspace(interval(1), interval(2), 300);
            
            [R, T] = meshgrid(r, thetas);

            func_vals = cellfun(@(c) c(R,T), allFunctions(used_indices), 'UniformOutput', false);

            vals = zeros(size(T));
            for iCoeff = 1:numel(func_vals)
                vals = vals + func_vals{iCoeff}*obj.Solution.Values(iCoeff);
            end

            X = R.*cos(T);
            Y = R.*sin(T);
            
            surf(X, Y, vals, "EdgeColor", "red")
        end

        function plotAirGap(obj, interval)
            arguments
                obj 
                interval = [0, 2*pi];
            end

            r = linspace(obj.CouplingInterface1.Radius, obj.CouplingInterface2.Radius, 10);
            thetas = linspace(interval(1), interval(2), 360);
            [R, T] = meshgrid(r, thetas);

            X = R.*cos(T);
            Y = R.*sin(T);
            vals = zeros(size(T));

            surf(X, Y, vals, "EdgeColor", "none", "FaceColor", TUDa_getColor("7b"))
        end
    end
end