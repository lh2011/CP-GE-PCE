classdef  EQS_System < EQS_Element
    properties
        Elements;
        TimeSolution;
        tSpan;
        NewtonMaxIter = 50;
        LinesearchMaxIter = 10;
        NewtonTolerance = 1e-12;
        Mlin, Klin, Flin;
    end

    methods
        function obj = EQS_System()
            obj.NumberDOF = 0;
            obj.resetElements();
        end

        function resetElements(obj)
            obj.Elements = struct('Element', {}, 'Indices', {}, 'IndicesExt', {}, 'IndicesInt', {});
        end

        function addElement(obj, element)
            pos = numel(obj.Elements) + 1;
            obj.Elements(pos).Element = element;

            obj.Elements(pos).IndicesInt = obj.NumberDOF + (1:element.NumberDOF);
            obj.Elements(pos).IndicesExt = obj.getExternalDOFs(element);
            obj.Elements(pos).Indices = [obj.Elements(pos).IndicesExt, obj.Elements(pos).IndicesInt];

            obj.NumberDOF = obj.NumberDOF + element.NumberDOF;
        end

        function [y0, Jac] = solveStatic(obj, t, y0, args)
            arguments
                obj 
                t = 0;
                y0 = zeros(obj.NumberDOF, 1);
                args.Verbose = true;
            end
            if numel(y0) ~= obj.NumberDOF
                error("Give an initial condition of the size of number of DOFs");
            end
            mytimer("Newton", "start");
            % Initialize matrices, save linear ones for efficiency
            obj.getBCmatrices();
            if isempty(obj.Klin)
                [~, obj.Klin] = obj.StiffnessMatrix(t, y0);
            end
            if isempty(obj.Flin)
                [~, obj.Flin] = obj.ForceVector(t, y0);
            end
            % Init
            K = obj.StiffnessMatrix(t, y0) + obj.Klin;
            f = obj.ForceVector(t, y0) + obj.Flin;
            % Apply boundary conditions
            b = obj.getBCvalues(t);
            Kred = K(obj.K, obj.K) + K(obj.K, obj.R)*obj.G + obj.G'*K(obj.R, obj.K) + obj.G'*K(obj.R, obj.R)*obj.G;
            fred = f(obj.K) + obj.G'*f(obj.R) + (K(obj.K, obj.R) + obj.G'*K(obj.R, obj.R))*obj.H*b;% - obj.G'*M(obj.R,obj.R)*obj.H*db;
            y0red = y0(obj.K);
            res = Kred*y0red + fred;
            rnorm = res'*res;
            for iNewton = 1:obj.NewtonMaxIter
                mytimer("Newton-Jacobian", "start");
                J = obj.JacobianMatrix(t, y0);
                Jac = Kred + J(obj.K, obj.K) + J(obj.K, obj.R)*obj.G + obj.G'*J(obj.R, obj.K) + obj.G'*J(obj.R, obj.R)*obj.G;
                mytimer("Newton-Jacobian", "stop");
                mytimer("Newton-Solve", "start")
                w = Jac\res;    %search direction
                mytimer("Newton-Solve", "stop")

                for iLinesearch = 0:obj.LinesearchMaxIter
                    tau = 0.5^iLinesearch;
                    y1red = y0red - tau*w;
                    % reconstruct
                    y1 = obj.reconstructSolution(t, y1red);
                    % Update system
                    mytimer("Newton-Res", "start");
                    K1 = obj.StiffnessMatrix(t, y1) + obj.Klin;
                    f1 = obj.ForceVector(t, y1) + obj.Flin;
                    % Apply boundary conditions
                    Kred1 = K1(obj.K, obj.K) + K1(obj.K, obj.R)*obj.G + obj.G'*K1(obj.R, obj.K) + obj.G'*K1(obj.R, obj.R)*obj.G;
                    fred1 = f1(obj.K) + obj.G'*f1(obj.R) + (K1(obj.K, obj.R) + obj.G'*K1(obj.R, obj.R))*obj.H*b;% - obj.G'*M(obj.R,obj.R)*obj.H*db;
                    % check residual
                    res1 = Kred1*y1red + fred1;
                    rnorm1 = res1'*res1;
                    mytimer("Newton-Res", "stop");
                    if rnorm1<rnorm || rnorm1 < obj.NewtonTolerance
                        Kred = Kred1;
                        res = res1;
                        rnorm = rnorm1;
                        y0 = y1;
                        y0red = y0(obj.K);
                        if args.Verbose
                            disp(['NewtonIter: ' 9 num2str(iNewton) 9 'Residual: ' 9 num2str(rnorm1) 9 'Linesearch ' num2str(iLinesearch)]);
                        end
                        break;
                    end
                end
                if iLinesearch == obj.LinesearchMaxIter
                    Kred = Kred1;
                    res = res1;
                    rnorm = rnorm1;
                    y0 = y1;
                    y0red = y0(obj.K);
                    if args.Verbose
                        disp(['NewtonIter: ' 9 num2str(iNewton) 9 'Residual: ' 9 num2str(rnorm1) 9 'Linesearch max value ' num2str(iLinesearch)]);
                    end
                end
                if rnorm1 < obj.NewtonTolerance
                    break
                end
            end
            y0 = obj.reconstructSolution(t, y0red);
            obj.postprocess(t, y0, zeros(size(y0)));
            mytimer("Newton", "stop");
        end

        function solveTransient(obj, tSpan, y0, solver, options)
            if ~exist("y0", "var") || isempty(y0)
                y0 = zeros(obj.NumberDOF, 1);
            end
            if numel(y0) ~= obj.NumberDOF
                error("Give an initial condition of the size of number of DOFs");
            end
            if ~exist("solver", "var") || isempty(solver)
                solver = @ode15s;
            end
            obj.getBCmatrices();
            y0 = y0(obj.K);

            M = obj.MassMatrixShrink(tSpan(1));

            standardOptions = odeset('InitialStep', 1e-3, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'Jacobian', @(t_, y_) obj.JacobianMatrixShrink(t_, y_), ...
                'Mass', M, 'MStateDependence', 'none', 'MassSingular', 'yes', 'Stats', 'on',  'OutputFcn', @(t_, y_, f_) obj.standardOutput(t_, y_, f_));
            standardOptions.("System") = obj;

            % Overwrite or extend the standard options
            if exist("options", "var")
                field_names = fieldnames (options);
                for iopts  = 1:numel (field_names)
                    if ~ isempty(options.(field_names{iopts}))
                        standardOptions.(field_names{iopts}) = options.(field_names{iopts});
                    end
                end
            end
            obj.tSpan = tSpan;
            obj.TimeSolution = solver(@(t, y) obj.solverFunction(t, y), tSpan, y0, standardOptions);
            % TBD self written solvers
            if strcmp(obj.TimeSolution.solver, 'odeIE') || strcmp(obj.TimeSolution.solver, 'odeTRA')
                yred = obj.TimeSolution.y;
                yredprime = gradient(obj.TimeSolution.y)./gradient(obj.TimeSolution.x);
            else
                [yred, yredprime] = deval(obj.TimeSolution, obj.TimeSolution.x);
            end

            [y, yprime] = obj.reconstructSolution(obj.TimeSolution.x, yred, yredprime);
            obj.postprocess(obj.TimeSolution.x, y, yprime);
        end

        function y0 = solveHarmonic(obj, t, y0, args)
            arguments
                obj 
                t = 0;
                y0 = zeros(obj.NumberDOF, 1);
                args.Verbose (1,1) logical = true;
                args.FixpointTolerance = 1e-6;
                args.PostprocessingSamples = 101;
                args.PostprocessingFrequencies = 1;
            end
            % Init
            W = obj.OmegaMatrix();
            % Apply boundary conditions
            obj.getBCmatrices();
            y0red = y0(obj.K);

            Wred = W(obj.K, obj.K);
            Kred = obj.StiffnessMatrixShrink(t, y0red);
            Mred = obj.MassMatrixShrink(t, y0red);
            Fred = obj.ForceVectorShrink(t, y0red);

            res = (-1i*Wred*Mred + Kred)*y0red + Fred;
            rnorm = res'*res;
            for iFixpoint = 1:50
                % Build new matrix and rhs
                y1red = (1i*Wred*Mred - Kred)\Fred;
                w = y1red - y0red;    %search direction

                for iLinesearch = 0:obj.LinesearchMaxIter
                    tau = 0.5^iLinesearch;
                    % update solution
                    y1red = y0red + tau*w;

                    % update system matrices
                    Kred1 = obj.StiffnessMatrixShrink(t, y1red);
                    Mred1 = obj.MassMatrixShrink(t, y1red);
                    Fred1 = obj.ForceVectorShrink(t, y1red);
                    % check residual
                    res1 = (-1i*Wred*Mred1 + Kred1)*y1red + Fred1;
                    rnorm1 = res1'*res1;
                    if rnorm1<rnorm || rnorm1 < args.FixpointTolerance
                        Kred = Kred1;
                        Mred = Mred1;
                        rnorm = rnorm1;
                        y0red = y1red;
                        if args.Verbose
                            disp(['FixpointIter: ' 9 num2str(iFixpoint) 9 'Residual: ' 9 num2str(rnorm1) 9 'Linesearch ' num2str(iLinesearch)]);
                        end
                        break;
                    end
                end
                if iLinesearch == obj.LinesearchMaxIter
                    Kred = Kred1;
                    Mred = Mred1;
                    rnorm = rnorm1;
                    y0red = y1red;
                    if args.Verbose
                        disp(['FixpointIter: ' 9 num2str(iFixpoint) 9 'Residual: ' 9 num2str(rnorm1) 9 'Linesearch max value ' num2str(iLinesearch)]);
                    end
                end
                if rnorm1 < args.FixpointTolerance
                    break
                end
            end

            % reconstruct for time solution
            omegas = diag(Wred);
            Tmax = max(2*pi./abs(omegas(omegas~=0)));
            t = linspace(0, Tmax*args.PostprocessingFrequencies, args.PostprocessingSamples);
            t(end) = [];
            yred = real(y0red.*exp(1i.*omegas.*t));
            yredprime = real(1i.*omegas.*y0red.*exp(1i.*omegas.*t));

            [y, yprime] = obj.reconstructSolution(t, yred, yredprime);
            obj.postprocess(t, y, yprime);
        end
        
        % Implementation according to P. Deuflhard: Newton Methods for Nonlinear Problems (2004) pp. 344-347
        % assumption: constant and time invariant mass matrix
        % T: time period length
        % nSample: Number of sample points for fourier estimation
        function y0red = solveHarmonicBalance(obj, T, y0red, harmonicsCos, harmonicsSin, nSample, verbose)
            if ~exist("T", "var")
                error("give the time interval for the fundamental!");
            end
            if ~exist("harmonicsCos", "var") || isempty(harmonicsCos)
                harmonicsCos = 1;
            end
            if ~exist("harmonicsSin", "var") || isempty(harmonicsSin)
                harmonicsSin = 1;
            end
            nCos = numel(harmonicsCos);
            nSin = numel(harmonicsSin);
            kCos = reshape(harmonicsCos, 1, []);
            kSin = reshape(harmonicsSin, 1, []);
            obj.getBCmatrices();    % initialize BC matrices;
            number_Kdofs = numel(obj.K);
            if ~exist("y0red", "var") || isempty(y0red)
                y0red = zeros(number_Kdofs, nCos+nSin);
            end
            if ~exist("nSample", "var") || isempty(nSample)
                nSample = 20;
            end
            if ~exist("verbose", "var") || isempty(verbose)
                verbose = true;
            end
            mytimer("HB", "start");
            % equidistant sampling points
            tp = linspace(0, T, nSample+1);
            tp = tp(1:end-1);
            omega = 2*pi/T;
            TwoT = 2/T;

            % Initialize matrices, save linear ones for efficiency
            if isempty(obj.Mlin)
                [~, obj.Mlin] = obj.MassMatrix(0, zeros(obj.NumberDOF, 1));
            end
            if isempty(obj.Klin)
                [~, obj.Klin] = obj.StiffnessMatrix(0, zeros(obj.NumberDOF, 1));
            end
            if isempty(obj.Flin)
                [~, obj.Flin] = obj.ForceVector(0, zeros(obj.NumberDOF, 1));
            end
            Mlinred = obj.Mlin(obj.K, obj.K) + obj.Mlin(obj.K, obj.R)*obj.G + obj.G'*obj.Mlin(obj.R, obj.K) + obj.G'*obj.Mlin(obj.R, obj.R)*obj.G;
            Klinred = obj.Klin(obj.K, obj.K) + obj.Klin(obj.K, obj.R)*obj.G + obj.G'*obj.Klin(obj.R, obj.K) + obj.G'*obj.Klin(obj.R, obj.R)*obj.G;
            Flinred = obj.Flin(obj.K) + obj.G'*obj.Flin(obj.R); 

            % indices for cos/sin coefficients
            ind_ak = reshape(1:number_Kdofs*nCos, number_Kdofs, nCos);
            ind_bk = reshape(number_Kdofs*nCos+(1:number_Kdofs*nSin), number_Kdofs, nSin);

            % discrete evaluations of fourier expansion
            % yred_tp = reshape(sum(y0red.*[cos(omega*kCos.*reshape(tp,1,1,[])), sin(omega*kSin.*reshape(tp,1,1,[]))], 2), numel(obj.K), numel(tp));
            yred_tp = reshape(sum(y0red.*[cospi(TwoT*kCos.*reshape(tp,1,1,[])), sinpi(TwoT*kSin.*reshape(tp,1,1,[]))], 2), numel(obj.K), numel(tp));
            y_tp = obj.reconstructSolution(tp, yred_tp);

            % evaluation of nonlinear state equation for sample points
            f_nonlin_tp = zeros(size(yred_tp));
            for itp = 1:nSample
                t = tp(itp);
                y = y_tp(:, itp);
                yred = yred_tp(:, itp);
                Knonlin = obj.StiffnessMatrix(t, y);
                Knonlinred = Knonlin(obj.K, obj.K) + Knonlin(obj.K, obj.R)*obj.G + obj.G'*Knonlin(obj.R, obj.K) + obj.G'*Knonlin(obj.R, obj.R)*obj.G;
                Fnonlin = obj.ForceVector(t, y);
                Fnonlinred = Fnonlin(obj.K) + obj.G'*Fnonlin(obj.R);
                % boundary conditions
                [b, db] = obj.getBCvalues(t);
                if any(b~=0, "all")
                    Kmat = obj.Klin + Knonlin;
                    Fnonlinred = Fnonlinred + (Kmat(obj.K, obj.R) + obj.G'*Kmat(obj.R, obj.R))*obj.H*b;
                end
                if any(db~=0, "all")
                    Fnonlinred = Fnonlinred - obj.G'*obj.Mlin(obj.R,obj.R)*obj.H*db;
                end
                f_nonlin_tp(:, itp) = Knonlinred*yred + Fnonlinred;
            end

            % fourier coefficients for cos terms
            alpha_k = zeros(number_Kdofs, nCos);
            for iCos = 1:nCos
                k = harmonicsCos(iCos);
                % linear contributions
                alpha_k(:, iCos) = Klinred*y0red(ind_ak(:, iCos));
                if k == 0
                    alpha_k(:, iCos) = 2*alpha_k(:, iCos) + 2*Flinred;  % times 2 for compatibility with 2/T for nonlinear
                end
                % nonlinear contributions
                alpha_k(:, iCos) = alpha_k(:, iCos) + 2/nSample*sum(f_nonlin_tp.*cos(k*omega*tp), 2);
            end
            % fourier coefficients for sin terms
            beta_k = zeros(number_Kdofs, nSin);
            for iSin = 1:nSin
                k = harmonicsSin(iSin);
                % linear contributions
                beta_k(:, iSin) = Klinred*y0red(ind_bk(:, iSin));
                % nonlinear contributions
                beta_k(:, iSin) = beta_k(:, iSin) + 2/nSample*sum(f_nonlin_tp.*sin(k*omega*tp), 2);
            end
            
            %%% generate residuum
            res = zeros(number_Kdofs*(nCos+nSin), 1);
            [~, commonIndCos, commonIndSin] = intersect(harmonicsCos, harmonicsSin);
            % contributions cos terms
            for iCos = 1:nCos
                k = harmonicsCos(iCos);
                if any(iCos == commonIndCos)
                    sinPosition = commonIndSin(find(iCos == commonIndCos, 1, "first"));
                    bk = y0red(ind_bk(:, sinPosition));
                else
                    bk = sparse(numel(obj.K), 1);
                end
                res(ind_ak(:, iCos)) = Mlinred*bk*omega*k - alpha_k(:, iCos);
            end
            % contributions sin terms
            for iSin = 1:nSin
                k = harmonicsSin(iSin);
                if any(iSin == commonIndSin)
                    cosPosition = commonIndCos(find(iSin == commonIndSin, 1, "first"));
                    ak = y0red(ind_ak(:, cosPosition));
                else
                    ak = sparse(numel(obj.K), 1);
                end
                res(ind_bk(:, iSin)) = Mlinred*ak*omega*k + beta_k(:, iSin);
            end

            rnorm = res'*res;

            for iNewton = 1:obj.NewtonMaxIter
                % Jacobian
                rows = cell(0, 1);
                cols = cell(0, 1);
                vals = cell(0, 1);
                % Jacobian (I): contributions from linear stiffness matrix
                [rwsKlin, clsKlin, valsKlin] = find(Klinred);
                for iCos = 1:nCos
                    rows{end+1,1} = ind_ak(rwsKlin, iCos);
                    cols{end+1,1} = ind_ak(clsKlin, iCos);
                    vals{end+1,1} = -valsKlin;
                    if harmonicsCos(iCos) == 0
                        vals{end,1} = vals{end,1}*2;
                    end
                end
                for iSin = 1:nSin
                    rows{end+1,1} = ind_bk(rwsKlin, iSin);
                    cols{end+1,1} = ind_bk(clsKlin, iSin);
                    vals{end+1,1} = valsKlin;
                end
                % Jacobian (II): contributions from nonlinearities
                for itp = 1:nSample
                    t = tp(itp);
                    y = y_tp(:, itp);
                    yred = yred_tp(:, itp);
                    J_nonlin = obj.StiffnessMatrix(t, y) + obj.JacobianMatrix(t, y);
                    Jred =  J_nonlin(obj.K, obj.K) + J_nonlin(obj.K, obj.R)*obj.G + obj.G'*J_nonlin(obj.R, obj.K) + obj.G'*J_nonlin(obj.R, obj.R)*obj.G;
                    [rwsJ, clsJ, valsJ] = find(Jred);
                    % for efficiency: check if nonlinearities occur at the same indices
                    if itp == 1
                        rwsJ1 = rwsJ;
                        clsJ1 = clsJ;
                        helpInd = numel(rows) + reshape(1:(nCos+nSin)^2, nCos+nSin, nCos+nSin)';
                    end
                    % for efficiency: check if nonlinearities occur at the same indices
                    if (itp>1) &&(numel(rwsJ) == numel(rwsJ1)) && (numel(clsJ) == numel(clsJ1)) && all(rwsJ==rwsJ1) && all(clsJ==clsJ1)
                        for iCos = 1:nCos
                            k1 = harmonicsCos(iCos);
                            for iCos2 = 1:nCos
                                k2 = harmonicsCos(iCos2);
                                vals{helpInd(iCos,iCos2)} = vals{helpInd(iCos,iCos2)} - 2/nSample*valsJ*cos(k1*omega*t)*cos(k2*omega*t);
                            end
                            for iSin2 = 1:nSin
                                k2 = harmonicsSin(iSin2);
                                vals{helpInd(iCos,nCos+iSin2)} = vals{helpInd(iCos,nCos+iSin2)} - 2/nSample*valsJ*cos(k1*omega*t)*sin(k2*omega*t);
                            end
                        end
                        for iSin = 1:nSin
                            k1 = harmonicsSin(iSin);
                            for iCos2 = 1:nCos
                                k2 = harmonicsCos(iCos2);
                                vals{helpInd(nCos+iSin,iCos2)} = vals{helpInd(nCos+iSin,iCos2)} + 2/nSample*valsJ*sin(k1*omega*t)*cos(k2*omega*t);
                            end
                            for iSin2 = 1:nSin
                                k2 = harmonicsSin(iSin2);
                                vals{helpInd(nCos+iSin,nCos+iSin2)} = vals{helpInd(nCos+iSin,nCos+iSin2)} + 2/nSample*valsJ*sin(k1*omega*t)*sin(k2*omega*t);
                            end
                        end
                    else    % otherwise generate new indices if it is different compared to before
                        for iCos = 1:nCos
                            k1 = harmonicsCos(iCos);
                            for iCos2 = 1:nCos
                                k2 = harmonicsCos(iCos2);
                                rows{end+1,1} = ind_ak(rwsJ, iCos);
                                cols{end+1,1} = ind_ak(clsJ, iCos2);
                                vals{end+1,1} = - 2/nSample*valsJ*cos(k1*omega*t)*cos(k2*omega*t);
                            end
                            for iSin2 = 1:nSin
                                k2 = harmonicsSin(iSin2);
                                rows{end+1,1} = ind_ak(rwsJ, iCos);
                                cols{end+1,1} = ind_bk(clsJ, iSin2);
                                vals{end+1,1} = - 2/nSample*valsJ*cos(k1*omega*t)*sin(k2*omega*t); 
                            end
                        end
                        for iSin = 1:nSin
                            k1 = harmonicsSin(iSin);
                            for iCos2 = 1:nCos
                                k2 = harmonicsCos(iCos2);
                                rows{end+1,1} = ind_bk(rwsJ, iSin);
                                cols{end+1,1} = ind_ak(clsJ, iCos2);
                                vals{end+1,1} = 2/nSample*valsJ*sin(k1*omega*t)*cos(k2*omega*t);
                            end
                            for iSin2 = 1:nSin
                                k2 = harmonicsSin(iSin2);
                                rows{end+1,1} = ind_bk(rwsJ, iSin);
                                cols{end+1,1} = ind_bk(clsJ, iSin2);
                                vals{end+1,1} = 2/nSample*valsJ*sin(k1*omega*t)*sin(k2*omega*t);
                            end
                        end
                    end
                end
                % Jacobian (III): contributions from mass matrix
                [rwsM, clsM, valsM] = find(Mlinred);
                for iHarm = 1:numel(commonIndCos)
                    k1 = harmonicsSin(commonIndSin(iHarm));
                    k2 = harmonicsCos(commonIndCos(iHarm));
                    assert(k1==k2);
                    iSin = commonIndSin(iHarm);
                    iCos = commonIndCos(iHarm);
                    rows{end+1,1} = ind_ak(rwsM,iCos);
                    cols{end+1,1} = ind_bk(clsM,iSin);
                    vals{end+1,1} = valsM*k1*omega;

                    rows{end+1,1} = ind_bk(rwsM,iSin);
                    cols{end+1,1} = ind_ak(clsM,iCos);
                    vals{end+1,1} = valsM*k1*omega;
                end
                mytimer("HB:Build", "start");
                Jac = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), number_Kdofs*(nCos+nSin), number_Kdofs*(nCos+nSin));
                mytimer("HB:Build", "stop");
                mytimer("HB:Solve", "start");
                w = Jac\res;    %search direction
                mytimer("HB:Solve", "stop");
                mytimer("HB:Linesearch", "start");
                for iLinesearch = 0:obj.LinesearchMaxIter
                    tau = 0.5^iLinesearch;
                    y1red = y0red - tau*reshape(w, numel(obj.K), nSin+nCos);

                    % discrete evaluations of fourier expansion
                    % yred_tp = reshape(sum(y1red.*[cos(omega*kCos.*reshape(tp,1,1,[])), sin(omega*kSin.*reshape(tp,1,1,[]))], 2), numel(obj.K), numel(tp));
                    yred_tp = reshape(sum(y1red.*[cospi(TwoT*kCos.*reshape(tp,1,1,[])), sinpi(TwoT*kSin.*reshape(tp,1,1,[]))], 2), numel(obj.K), numel(tp));
                    y_tp = obj.reconstructSolution(tp, yred_tp);
        
                    % evaluation of nonlinear state equation for sample points
                    f_nonlin_tp = zeros(size(yred_tp));
                    for itp = 1:nSample
                        t = tp(itp);
                        y = y_tp(:, itp);
                        yred = yred_tp(:, itp);
                        mytimer("HB:StiffnessMatrix", "start")
                        Knonlin = obj.StiffnessMatrix(t, y);
                        mytimer("HB:StiffnessMatrix", "stop")
                        Knonlinred = Knonlin(obj.K, obj.K) + Knonlin(obj.K, obj.R)*obj.G + obj.G'*Knonlin(obj.R, obj.K) + obj.G'*Knonlin(obj.R, obj.R)*obj.G;
                        Fnonlin = obj.ForceVector(t, y);
                        Fnonlinred = Fnonlin(obj.K) + obj.G'*Fnonlin(obj.R);
                        % boundary conditions
                        [b, db] = obj.getBCvalues(t);
                        if any(b~=0, "all")
                            Kmat = obj.Klin + Knonlin;
                            Fnonlinred = Fnonlinred + (Kmat(obj.K, obj.R) + obj.G'*Kmat(obj.R, obj.R))*obj.H*b;
                        end
                        if any(db~=0, "all")
                            Fnonlinred = Fnonlinred - obj.G'*obj.Mlin(obj.R,obj.R)*obj.H*db;
                        end
                        f_nonlin_tp(:, itp) = Knonlinred*yred + Fnonlinred;
                    end
        
                    % fourier coefficients for cos terms
                    alpha_k = zeros(number_Kdofs, nCos);
                    for iCos = 1:nCos
                        k = harmonicsCos(iCos);
                        % linear contributions
                        alpha_k(:, iCos) = Klinred*y1red(ind_ak(:, iCos));
                        if k == 0
                            alpha_k(:, iCos) = 2*alpha_k(:, iCos) + 2*Flinred;  % times 2 for compatibility with 2/T for nonlinear
                        end
                        % nonlinear contributions
                        alpha_k(:, iCos) = alpha_k(:, iCos) + 2/nSample*sum(f_nonlin_tp.*cos(k*omega*tp), 2);
                    end
                    % fourier coefficients for sin terms
                    beta_k = zeros(number_Kdofs, nSin);
                    for iSin = 1:nSin
                        k = harmonicsSin(iSin);
                        % linear contributions
                        beta_k(:, iSin) = Klinred*y1red(ind_bk(:, iSin));
                        % nonlinear contributions
                        beta_k(:, iSin) = beta_k(:, iSin) + 2/nSample*sum(f_nonlin_tp.*sin(k*omega*tp), 2);
                    end
                    
                    %%% generate residuum
                    res1 = zeros(number_Kdofs*(nCos+nSin), 1);
                    [~, commonIndCos, commonIndSin] = intersect(harmonicsCos, harmonicsSin);
                    % contributions cos terms
                    for iCos = 1:nCos
                        k = harmonicsCos(iCos);
                        if any(iCos == commonIndCos)
                            sinPosition = commonIndSin(find(iCos == commonIndCos, 1, "first"));
                            bk = y1red(ind_bk(:, sinPosition));
                        else
                            bk = sparse(numel(obj.K), 1);
                        end
                        res1(ind_ak(:, iCos)) = Mlinred*bk*omega*k - alpha_k(:, iCos);
                    end
                    % contributions sin terms
                    for iSin = 1:nSin
                        k = harmonicsSin(iSin);
                        if any(iSin == commonIndSin)
                            cosPosition = commonIndCos(find(iSin == commonIndSin, 1, "first"));
                            ak = y1red(ind_ak(:, cosPosition));
                        else
                            ak = sparse(numel(obj.K), 1);
                        end
                        res1(ind_bk(:, iSin)) = Mlinred*ak*omega*k + beta_k(:, iSin);
                    end
                    rnorm1 = res1'*res1;

                    if rnorm1<rnorm || rnorm1 < obj.NewtonTolerance
                        res = res1;
                        rnorm = rnorm1;
                        y0red = y1red;
                        if verbose
                            disp(['NewtonIter: ' 9 num2str(iNewton) 9 'Residual: ' 9 num2str(rnorm1) 9 'Linesearch ' num2str(iLinesearch)]);
                        end
                        break;
                    end
                end
                mytimer("HB:Linesearch", "stop");
                if iLinesearch == obj.LinesearchMaxIter
                    res = res1;
                    rnorm = rnorm1;
                    y0red = y1red;
                    if verbose
                        disp(['NewtonIter: ' 9 num2str(iNewton) 9 'Residual: ' 9 num2str(rnorm1) 9 'Linesearch max value ' num2str(iLinesearch)]);
                    end
                end
                if rnorm1 < obj.NewtonTolerance
                    break
                end 
            end

            yred_tp = reshape(sum(y0red.*[cos(omega*kCos.*reshape(tp,1,1,[])), sin(omega*kSin.*reshape(tp,1,1,[]))], 2), numel(obj.K), numel(tp));
            yredprime_tp = reshape(sum(y0red.*[-kCos.*sin(omega*kCos.*reshape(tp,1,1,[])), kSin.*cos(omega*kSin.*reshape(tp,1,1,[]))], 2), numel(obj.K), numel(tp));
            
            [ytp, yprimetp] = obj.reconstructSolution(tp, yred_tp, yredprime_tp);
            obj.postprocess(tp, ytp, yprimetp);
            mytimer("HB", "stop");
        end
        
        function stop = standardOutput(obj, t, y, flag)
            stop = false;
            switch flag
                case 'init'
                    progressbar('Starting Integration: ')
                case 'done'
                    progressbar('finished')
                otherwise
                    progressbar((t-obj.tSpan(1))/(obj.tSpan(end)-obj.tSpan(1)));
            end
        end

        function res = solverFunction(obj, t, yred)
            Kred = obj.StiffnessMatrixShrink(t, yred);
            Fred = obj.ForceVectorShrink(t, yred);

            res = Kred*yred + Fred;
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            if nargout == 1
                % Mnonlin = sparse(obj.NumberDOF, obj.NumberDOF);
                rows = cell(numel(obj.Elements), 1);
                cols = cell(numel(obj.Elements), 1);
                vals = cell(numel(obj.Elements), 1);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    Mnonlin_iel = obj.Elements(iEl).Element.MassMatrix(t, y_iel);
                    % Mnonlin(ind_iel, ind_iel) = Mnonlin(ind_iel, ind_iel) + Mnonlin_iel;
                    [r, c, v] = find(Mnonlin_iel);
                    rows{iEl} = reshape(ind_iel(r), [], 1);
                    cols{iEl} = reshape(ind_iel(c), [], 1);
                    vals{iEl} = reshape(v, [], 1);
                end
                Mnonlin = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), obj.NumberDOF, obj.NumberDOF);
            else
                % Mnonlin = sparse(obj.NumberDOF, obj.NumberDOF);
                % Mlin = sparse(obj.NumberDOF, obj.NumberDOF);
                rowsLin = cell(numel(obj.Elements), 1);
                colsLin = cell(numel(obj.Elements), 1);
                valsLin = cell(numel(obj.Elements), 1);
                rowsNonlin = cell(numel(obj.Elements), 1);
                colsNonlin = cell(numel(obj.Elements), 1);
                valsNonlin = cell(numel(obj.Elements), 1);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    [Mnonlin_iel, Mlin_iel] = obj.Elements(iEl).Element.MassMatrix(t, y_iel);
                    % Mnonlin(ind_iel, ind_iel) = Mnonlin(ind_iel, ind_iel) + Mnonlin_iel;
                    % Mlin(ind_iel, ind_iel) = Mlin(ind_iel, ind_iel) + Mlin_iel;
                    [rNL, cNL, vNL] = find(Mnonlin_iel);
                    rowsNonlin{iEl} = reshape(ind_iel(rNL), [], 1);
                    colsNonlin{iEl} = reshape(ind_iel(cNL), [], 1);
                    valsNonlin{iEl} = reshape(vNL, [], 1);
                    [rL, cL, vL] = find(Mlin_iel);
                    rowsLin{iEl} = reshape(ind_iel(rL), [], 1);
                    colsLin{iEl} = reshape(ind_iel(cL), [], 1);
                    valsLin{iEl} = reshape(vL, [], 1);
                end
                Mlin = sparse(cell2mat(rowsLin), cell2mat(colsLin), cell2mat(valsLin), obj.NumberDOF, obj.NumberDOF);
                Mnonlin = sparse(cell2mat(rowsNonlin), cell2mat(colsNonlin), cell2mat(valsNonlin), obj.NumberDOF, obj.NumberDOF);
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            if nargout == 1
                % Knonlin = sparse(obj.NumberDOF, obj.NumberDOF);
                rows = cell(numel(obj.Elements), 1);
                cols = cell(numel(obj.Elements), 1);
                vals = cell(numel(obj.Elements), 1);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    Knonlin_iel = obj.Elements(iEl).Element.StiffnessMatrix(t, y_iel);
                    % Knonlin(ind_iel, ind_iel) = Knonlin(ind_iel, ind_iel) + Knonlin_iel;
                    [r, c, v] = find(Knonlin_iel);
                    rows{iEl} = reshape(ind_iel(r), [], 1);
                    cols{iEl} = reshape(ind_iel(c), [], 1);
                    vals{iEl} = reshape(v, [], 1);
                end
                Knonlin = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), obj.NumberDOF, obj.NumberDOF);
            else
                rowsLin = cell(numel(obj.Elements), 1);
                colsLin = cell(numel(obj.Elements), 1);
                valsLin = cell(numel(obj.Elements), 1);
                rowsNonlin = cell(numel(obj.Elements), 1);
                colsNonlin = cell(numel(obj.Elements), 1);
                valsNonlin = cell(numel(obj.Elements), 1);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    [Knonlin_iel, Klin_iel] = obj.Elements(iEl).Element.StiffnessMatrix(t, y_iel);
                    % Knonlin(ind_iel, ind_iel) = Knonlin(ind_iel, ind_iel) + Knonlin_iel;
                    % Klin(ind_iel, ind_iel) = Klin(ind_iel, ind_iel) + Klin_iel;
                    [rNL, cNL, vNL] = find(Knonlin_iel);
                    rowsNonlin{iEl} = reshape(ind_iel(rNL), [], 1);
                    colsNonlin{iEl} = reshape(ind_iel(cNL), [], 1);
                    valsNonlin{iEl} = reshape(vNL, [], 1);
                    [rL, cL, vL] = find(Klin_iel);
                    rowsLin{iEl} = reshape(ind_iel(rL), [], 1);
                    colsLin{iEl} = reshape(ind_iel(cL), [], 1);
                    valsLin{iEl} = reshape(vL, [], 1);
                end
                Klin = sparse(cell2mat(rowsLin), cell2mat(colsLin), cell2mat(valsLin), obj.NumberDOF, obj.NumberDOF);
                Knonlin = sparse(cell2mat(rowsNonlin), cell2mat(colsNonlin), cell2mat(valsNonlin), obj.NumberDOF, obj.NumberDOF);
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            if nargout == 1
                Fnonlin = zeros(obj.NumberDOF, 1);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    Fnonlin_iel = obj.Elements(iEl).Element.ForceVector(t, y_iel);
                    Fnonlin(ind_iel) = Fnonlin(ind_iel) + Fnonlin_iel;
                end
            else
                Fnonlin = zeros(obj.NumberDOF, 1);
                Flin = zeros(obj.NumberDOF, 1);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    [Fnonlin_iel, Flin_iel] = obj.Elements(iEl).Element.ForceVector(t, y_iel);
                    Fnonlin(ind_iel) = Fnonlin(ind_iel) + Fnonlin_iel;
                    Flin(ind_iel) = Flin(ind_iel) + Flin_iel;
                end
            end
        end

        function J = JacobianMatrix(obj, t, y)
            % J = sparse(obj.NumberDOF, obj.NumberDOF);
            rows = cell(numel(obj.Elements), 1);
            cols = cell(numel(obj.Elements), 1);
            vals = cell(numel(obj.Elements), 1);
            for iEl = 1:numel(obj.Elements)
                ind_iel = obj.Elements(iEl).Indices;
                y_iel = y(ind_iel);
                J_iel = obj.Elements(iEl).Element.JacobianMatrix(t, y_iel);
                [r, c, v] = find(J_iel);
                % J(ind_iel, ind_iel) = J(ind_iel, ind_iel) + J_iel;
                rows{iEl} = reshape(ind_iel(r), [], 1);
                cols{iEl} = reshape(ind_iel(c), [], 1);
                vals{iEl} = reshape(v, [], 1);
            end
            J = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), obj.NumberDOF, obj.NumberDOF);
        end

        % set omega for all elements for harmonic solution
        function setOmega(obj, omega)
            obj.omega = omega;
            for iEl = 1:numel(obj.Elements)
                obj.Elements(iEl).Element.setOmega(omega);
            end
        end

        % assemble diagonal matrix for all omega values
        function W = OmegaMatrix(obj)
            W = sparse(obj.NumberDOF, obj.NumberDOF);
            for iEl = 1:numel(obj.Elements)
                ind_iel = obj.Elements(iEl).IndicesInt;
                W_iel = obj.Elements(iEl).Element.OmegaMatrix();
                W(ind_iel, ind_iel) = W(ind_iel, ind_iel) + W_iel;
            end
        end

        function [G, H] = getBCmatrices(obj)
            % values were already called and initialized
            if issparse(obj.G) && issparse(obj.H)
                G = obj.G;
                H = obj.H;
                return
            end
            G = sparse(obj.NumberDOF, obj.NumberDOF);    % G = (RxK)
            H = sparse(obj.NumberDOF, obj.NumberDOF);    % H = (RxR)
            for iEl = 1:numel(obj.Elements)
                checkEl = obj.Elements(iEl).Element;
                [G_iel, H_iel] = checkEl.getBCmatrices();
                [localDofsR, localDofsK] = checkEl.getBCindices();
                globalDofsR = obj.Elements(iEl).IndicesInt(localDofsR);
                globalDofsK = obj.Elements(iEl).IndicesInt(localDofsK);
                G(globalDofsR, globalDofsK) = G(globalDofsR, globalDofsK) + G_iel;
                H(globalDofsR, globalDofsR) = H(globalDofsR, globalDofsR) + H_iel;
            end
            [dofsR, dofsK] = obj.getBCindices();
            G = G(dofsR, dofsK);
            H = H(dofsR, dofsR);
            obj.G = G;
            obj.H = H;
        end

        function [dofsR, dofsK] = getBCindices(obj)
            % values were already called and initialized
            if ~isempty(obj.K) || ~isempty(obj.R)
                dofsR = obj.R;
                dofsK = obj.K;
                return
            end
            % otherwise: go through elements and collect BC indices
            dofsR = [];
            dofsK = [];
            for iEl = 1:numel(obj.Elements)
                checkEl = obj.Elements(iEl).Element;
                globalDofs = obj.Elements(iEl).IndicesInt;
                [localDofsR, localDofsK] = checkEl.getBCindices();
                globalDofsR = globalDofs(localDofsR);
                globalDofsK = globalDofs(localDofsK);
                dofsR = [dofsR, globalDofsR];
                dofsK = [dofsK, globalDofsK];
            end
            if isempty(dofsR)
                dofsR = zeros(0,1); % fix if there are no BC
            end
            obj.R = dofsR;
            obj.K = dofsK;
            assert(numel(dofsR) == numel(unique(dofsR)))
            assert(numel(dofsK) == numel(unique(dofsK)))
        end

        function [b, db] = getBCvalues(obj, t)
            b = sparse(obj.NumberDOF, 1);
            db = sparse(obj.NumberDOF, 1);
            for iEl = 1:numel(obj.Elements)
                checkEl = obj.Elements(iEl).Element;
                [b_iel, db_iel] = checkEl.getBCvalues(t);
                localDofsR = checkEl.getBCindices();
                globalDofsR = obj.Elements(iEl).IndicesInt(localDofsR);
                b(globalDofsR) = b(globalDofsR) + b_iel;
                db(globalDofsR) = db(globalDofsR) + db_iel;
            end
            dofsR = obj.getBCindices();
            b = b(dofsR, 1);
            db = db(dofsR, 1);
        end

        function extDofs = getExternalDOFs(obj, element)
            extDofs = [];
            for iElExt = 1:numel(element.ExternalElements)
                extEl = element.ExternalElements(iElExt).Element;
                extInd = element.ExternalElements(iElExt).LocalIndices;
                extDofs_iel = obj.getIndices(extEl, extInd);
                if isscalar(extDofs_iel) && extDofs_iel == -1
                    error('EQS_System:ElementNotYetAdded', 'The external element %s is not found. Add it first to the system.', element.Name);
                end
                extDofs = [extDofs, extDofs_iel];
            end
            assert(all(extDofs == unique(extDofs, 'stable')));
        end

        function elementIndices = getIndices(obj, element, localIndices)
            if ~exist("localIndices", "var")
                localIndices = 1:element.NumberDOF;
            end
            elementIndices = -1;
            for iEl = 1:numel(obj.Elements)
                checkEl = obj.Elements(iEl).Element;
                globalDofs = obj.Elements(iEl).IndicesInt;
                if isa(checkEl, "EQS_System")
                    if checkEl == element
                        elementIndices = globalDofs(localIndices);
                    else
                        systemIndices = checkEl.getIndices(element, localIndices);
                        if ~isscalar(systemIndices) || systemIndices~=-1
                            elementIndices = globalDofs(systemIndices);
                        end
                    end
                elseif isa(checkEl, "EQS_Element")
                    if checkEl == element
                        elementIndices = globalDofs(localIndices);
                    end
                else
                    warning("You need either a EQS_Element or EQS_System here");
                end
                if elementIndices ~= -1
                    break;
                end
            end
        end

        function redIndices = getIndicesReduced(obj, element, localIndices)
            if ~exist("localIndices", "var")
                localIndices = 1:element.NumberDOF;
            end
            [~, K] = obj.getBCindices();
            elementIndices = obj.getIndices(element, localIndices);
            [~, redIndices] = intersect(K, elementIndices);
        end

        function [maxErr, Jac_num, Jac_ana] = checkJacobian(obj, t, y0red)
            if ~exist("t", "var")
                t = 0;
            end
            [dofsR, dofsK] = obj.getBCindices();
            if ~exist("y0red", "var")
                y0red = zeros(numel(dofsK), 1);
            end
            Step = 1e-8;
            f0 = obj.solverFunction(t, y0red);
            Jac_ana = full(obj.JacobianMatrixShrink(t, y0red));
            Jac_num = zeros(size(Jac_ana));
            progressbar('testing Jac')
            for idof = 1:numel(dofsK)
                y1 = y0red;
                y1(idof) = y1(idof) + Step;
                f1 = obj.solverFunction(t, y1);
                Jac_num(:, idof) = (f1-f0)/Step;
                progressbar(idof/numel(dofsK));
            end
            progressbar('done')
            e_rel = 1e-2;
            e_abs = 1e-3;
            err = abs(Jac_num-Jac_ana)./(max(abs(Jac_ana*e_rel), e_abs*ones(size(Jac_ana))));

            [rows, cols] = find(err>1);
            ind = sub2ind(size(Jac_ana), rows, cols);
            a = Jac_ana(ind);
            b = Jac_num(ind);
            maxErr = max(err, [], 'all');
        end

        % function plotElements(obj)
        %     figure(123)
        %     elements = getElementList(obj);
        %     mat = zeros(obj.NumberDOF, obj.NumberDOF, 'uint32');
        %     for iEl = 1:numel(elements)
        %         mat(elements(iEl).indicesInt, elements(iEl).indicesInt) = iEl;
        %         mat(elements(iEl).indicesInt, elements(iEl).indicesExt) = iEl;
        %         mat(elements(iEl).indicesExt, elements(iEl).indicesInt) = iEl;
        %     end
        %     imagesc(mat)
        % end
        %
        % function elements = getElementList(obj)
        %     elements = struct('element', {}, 'indicesExt', {}, 'indicesInt', {});
        %     for iEl = 1:numel(obj.Elements)
        %         elementsEl = obj.Elements(iEl).Element.getElementList();
        %         for iEl2 = 1:numel(elementsEl)
        %             elements(end+1).element = elementsEl(iEl2).element;
        %             elements(end).indicesInt = obj.Elements(iEl).Indices(elementsEl(iEl2).indicesInt);
        %             elements(end).indicesExt = obj.Elements(iEl).Indices(elementsEl(iEl2).indicesExt);
        %         end
        %     end
        % end

        %%%%%%%% postprocessing

        function postprocess(obj, ts, ys, ysprime)
            obj.Solution.Time = ts;
            obj.Solution.Values = ys;
            if ~exist("ysprime", "var") || isempty(ysprime)
                ysprime = gradient(ys)./gradient(reshape(ts, 1, []));
            end
            obj.Solution.ValuesDer = ysprime;
            for iEl = 1:numel(obj.Elements)
                obj.Elements(iEl).Element.postprocess(ts, ys(obj.Elements(iEl).Indices, :), ysprime(obj.Elements(iEl).Indices, :));
            end
        end
    end
end
