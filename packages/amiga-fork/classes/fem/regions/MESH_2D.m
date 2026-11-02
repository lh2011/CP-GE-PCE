% generic mesh for integrator regimes

classdef MESH_2D

    properties
        rdim;
        ndim;
        nqn;
        nel;
        geo_map;
        geo_map_jac;
        quad_weights;
        quad_nodes;
        jacdet;
        length;
        map, map_jac, map_det;
    end

     methods
        function obj = MESH_2D(meshEval, args)
            arguments
                meshEval,
                args.length = 1;
            end
            jacdet_cell = cellfun(@(s) s.jacdet, meshEval, 'UniformOutput', false);
            quad_nodes_cell = cellfun(@(s) s.quad_nodes, meshEval, 'UniformOutput', false);
            quad_weights_cell = cellfun(@(s) s.quad_weights, meshEval, 'UniformOutput', false);
            geo_map_cell = cellfun(@(s) s.geo_map, meshEval, 'UniformOutput', false);
            geo_map_jac_cell = cellfun(@(s) s.geo_map_jac, meshEval, 'UniformOutput', false);
            nel_cell = cellfun(@(s) s.nel, meshEval, 'UniformOutput', false);

            % condense meshes
            obj.jacdet = cat(2, jacdet_cell{:});
            obj.quad_nodes = cat(3, quad_nodes_cell{:});
            obj.quad_weights = cat(2, quad_weights_cell{:});
            obj.geo_map = cat(3, geo_map_cell{:});
            obj.geo_map_jac = cat(4, geo_map_jac_cell{:});
            obj.length = args.length;
            obj.nel = sum(cat(1, nel_cell{:}));
            obj.nqn = meshEval{1}.nqn; % assume that all elements have same nqn
            obj.rdim = 2; % assume 2D
            obj.ndim = 2; % assume 2D
        end
     end
end
