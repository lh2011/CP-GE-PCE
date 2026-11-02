% generic space for integrator regimes

classdef SPACE_2D
    properties
        connectivity;
        ndof;
        % dofs_ornt = [];
        gnum
        nsh_max;
        ncomp;
        rdim = 2;
        bsx_func;
        bsx_func_der;
        shape_functions
        shape_function_gradients;
    end

    methods
        function obj = SPACE_2D(spaceEval, spacesConnectivity, ncomp)
            obj.ncomp = ncomp;
            shape_functions_cell = cellfun(@(s) s.shape_functions, spaceEval, 'UniformOutput', false);
            shape_function_gradients_cell = cellfun(@(s) s.shape_function_gradients, spaceEval, 'UniformOutput', false);
            connectivity_cell = cellfun(@(s) s.connectivity, spaceEval, 'UniformOutput', false);
            for iConn = 1:numel(spacesConnectivity)
                % convert local to global connectivity
                connectivity_cell{iConn} = spacesConnectivity{iConn}(connectivity_cell{iConn});
            end

            % depending on ncomp, the indices are different
            if ncomp == 1
                indNqn = 1;
                indNsh = 2;
                indShp = 3;
                indShpGrad = 4;
            else
                indNqn = 2;
                indNsh = 3;
                indShp = 4;
                indShpGrad = 5;
            end

            maxNqn = max(cellfun(@(c) size(c, indNqn), shape_functions_cell));
            maxNsh = max(cellfun(@(c) size(c, indNsh), shape_functions_cell));
            % quick fix in case there is a different number of  shape functions
            % (e.g. for NURBS geometry functions): fill missing with zeros
            for iSpaces = 1:numel(shape_functions_cell)
                sizeNqn = size(shape_functions_cell{iSpaces}, indNqn);
                sizeNsh = size(shape_functions_cell{iSpaces}, indNsh);
                if sizeNqn ~= maxNqn
                    error("SPACE_2D: Only equal number of quadrature points per element considered!");
                end
                if sizeNsh < maxNsh
                    connectivity_cell{iSpaces}(sizeNsh+1:maxNsh,:) = connectivity_cell{iSpaces}(1,1);
                    if ncomp == 1
                        shape_functions_cell{iSpaces}(:,sizeNsh+1:maxNsh,:) = 0;
                        shape_function_gradients_cell{iSpaces}(:,:,sizeNsh+1:maxNsh,:) = 0;
                    else
                        shape_functions_cell{iSpaces}(:,:,sizeNsh+1:maxNsh,:) = 0;
                        shape_function_gradients_cell{iSpaces}(:,:,:,sizeNsh+1:maxNsh,:) = 0;
                    end
                end
            end

            obj.shape_functions = cat(indShp, shape_functions_cell{:});
            obj.shape_function_gradients = cat(indShpGrad, shape_function_gradients_cell{:});
            globalConn = cat(2, connectivity_cell{:});
            % reduce global connectivity to only DOFs that are actually used
            obj.connectivity = zeros(size(globalConn));
            [obj.gnum, ~, obj.connectivity(:)] = unique(globalConn);
            obj.ndof = numel(obj.gnum);
            obj.nsh_max = size(obj.shape_functions, indNsh);
        end
    end
end
