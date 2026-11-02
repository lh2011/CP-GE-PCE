classdef MESH_2Dtri

    properties
        rdim;
        nqn;
        nel;
        geo_map;
        geo_map_jac;
        quad_weights;
        quad_nodes;
        jacdet;
        length;
        map, map_jac, map_det;
        map_inv;

        quad_param
    end

    methods
        function obj = MESH_2Dtri(nodes, elements, args)
            arguments
                nodes (:,2) double;
                elements (:,3) double;
                args.DegreeQuadrature {mustBeMember(args.DegreeQuadrature, [0,1,2])} = 1;
                args.length = 1;
                args.qn = [];   % manual setting
                args.qw = [];   % manual setting
            end

            switch args.DegreeQuadrature
                case 0
                    nqn = 1;
                    qn = {1/3, 1/3};
                    qw = 1/2;
                case 1
                    nqn = 3;
                    qn = {[1/6, 2/3, 1/6], [1/6, 1/6, 2/3]};
                    qw = [1/6; 1/6; 1/6];
                case 2
                    nqn = 4;
                    qn = {[1/3, 1/5, 3/5, 1/5], [1/3, 1/5, 1/5, 3/5]};
                    qw = [-27/96; 25/96; 25/96; 25/96];
                otherwise
                    error('MESH_2Dtri:Quadrature rule not yet implemented with degree %i', args.DegreeQuadrature);
            end
            if ~isempty(args.qn) && ~isempty(args.qw)
                qn = args.qn;
                qw = args.qw;
                nqn = numel(qw);
            end

            obj.length = args.length;

            obj.rdim = 2;
            obj.nqn = nqn;
            obj.nel = size(elements, 1);

            ax = nodes(elements(:, 1), 1);
            ay = nodes(elements(:, 1), 2);
            bx = nodes(elements(:, 2), 1);
            by = nodes(elements(:, 2), 2);
            cx = nodes(elements(:, 3), 1);
            cy = nodes(elements(:, 3), 2);
            vals = permute(cat(3, [ax-cx, bx-cx], [ay-cy, by-cy]), [3,2,1]);
            rhs = permute([cx, cy], [2,1]);
            % PTS: 2xnpts (x nel)
            obj.map = @(PTS) permute(sum(vals.*permute(PTS, [4,1,3,2]), 2), [1,4,3,2]) + permute(rhs, [1,3,2]);
            obj.map_jac = @(PTS) repmat(permute(vals, [1,2,4,3]), [1, 1, size(PTS, 2)]);
            obj.map_det = @(PTS) repmat(permute(abs((ax-cx).*(by-cy) - (ay-cy).*(bx-cx)), [2,1]), [size(PTS, 2), 1]);
            
            % obj.map_inv{iEl} = @(PTS) [by-cy, -(bx-cx); -(ay-cy), ax-cx]/((ax-cx)*(by-cy) - (ay-cy)*(bx-cx)) *(PTS - [cx;cy]);
            lhs = permute(cat(3, [by-cy, cy-ay], [cx-bx, ax-cx]) ./ ((ax-cx).*(by-cy) - (ay-cy).*(bx-cx)), [2,3,1]);
            obj.map_inv = @(PTS) permute(sum(lhs.*permute((PTS - permute(cat(3, cx,cy), [3,2,1])), [4,1,3,2]), 2), [1,4,3,2]);

            obj.quad_param = [qn{1}; qn{2}];
            obj.geo_map = obj.map(obj.quad_param);
            obj.geo_map_jac = obj.map_jac(obj.quad_param);
            obj.quad_weights = repmat(qw, [1, obj.nel]);
            obj.quad_nodes = repmat(obj.quad_param, [1, 1, obj.nel]);
            obj.jacdet = obj.map_det(obj.quad_param);
        end
    end
end
