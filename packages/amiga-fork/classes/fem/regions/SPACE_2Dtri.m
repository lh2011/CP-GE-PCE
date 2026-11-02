% space for integrator regimes

classdef SPACE_2Dtri

    properties
        connectivity;
        ndof;
        % dofs_ornt = [];
        gnum
        nsh_max;
        % degree = 1;
        ncomp = 1; % TBD get this when importing
        rdim = 2;
        bsx_func;
        bsx_func_der;
        shape_functions
        shape_function_gradients;
        degree;
    end

    methods
        function obj = SPACE_2Dtri(mesh2D, elements, args)
            arguments
                mesh2D MESH_2Dtri;
                elements (:,3) double;
                args.degree {mustBeMember(args.degree, 1)} = 1; % only degree 1 for now
            end

            obj.connectivity = zeros(3, mesh2D.nel);
            [obj.gnum, ~,  obj.connectivity(:)] = unique(elements');
            obj.ndof = numel(obj.gnum);
            obj.nsh_max = 3;
            obj.degree = args.degree;
            obj.ncomp = 1;

            obj.bsx_func = @(xi, eta) cat(2, xi, eta, 1 - xi - eta);
            obj.bsx_func_der = @(xi, eta) cat(1, cat(3, ones(size(xi)), zeros(size(xi)), -ones(size(xi))),...
                cat(3, zeros(size(xi)), ones(size(xi)), -ones(size(xi))));
            xi = permute(mesh2D.quad_nodes(1,:,:), [2,1,3]);
            eta = permute(mesh2D.quad_nodes(2,:,:), [2,1,3]);
            obj.shape_functions = obj.bsx_func(xi, eta);
            % gradNihat: dim x nqn x N x nel
            gradNihat = obj.bsx_func_der(permute(xi, [2,1,4,3]), permute(eta, [2,1,4,3]));
            % Jinv: dim x dim x nqn x nel
            Jinv = pageinv(mesh2D.map_jac(mesh2D.quad_param));
            obj.shape_function_gradients = sum(permute(Jinv, [2,3,5,4,1]).*permute(gradNihat, [5,2,3,4,1]), 5);
        end
    end
end
