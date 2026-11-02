classdef BC_PeriodicNormal < BC
    properties (Access = public)
        LeftDOFs, RightDOFs;
        LeftDOFsX, LeftDOFsY;
        RightDOFsX, RightDOFsY;
    end

    methods (Access = public)
        function obj = BC_PeriodicNormal(parentElement, leftDOFsX, leftDOFsY, rightDOFsX, rightDOFsY, nL, nR)
            leftDOFsX = reshape(leftDOFsX, 1, []); leftDOFsY = reshape(leftDOFsY, 1, []);
            rightDOFsX = reshape(rightDOFsX, 1, []); rightDOFsY = reshape(rightDOFsY, 1, []);
            assert(isempty(intersect(union(leftDOFsX,leftDOFsX), union(leftDOFsY,leftDOFsY))), "Left and right DOFs must not be identical!");
            % Reduce DOFs for possible multiplicity with other Boundary Conditions
            obj@BC(parentElement, union([leftDOFsX, leftDOFsY], [rightDOFsX, rightDOFsY]));
            obj.LeftDOFsX = intersect(leftDOFsX, obj.BoundaryDOFs, 'stable');
            obj.LeftDOFsY = intersect(leftDOFsY, obj.BoundaryDOFs, 'stable');
            obj.RightDOFsX = intersect(rightDOFsX, obj.BoundaryDOFs, 'stable');
            obj.RightDOFsY = intersect(rightDOFsY, obj.BoundaryDOFs, 'stable');
            obj.LeftDOFs = [obj.LeftDOFsX, obj.LeftDOFsY];
            obj.RightDOFs = [obj.RightDOFsX, obj.RightDOFsY];

            % only set LeftDOFsX as the remove DOFs
            if nL(1) ~= 0
                obj.RemoveDOFs = obj.LeftDOFsX;
            elseif nL(2) ~= 0
                obj.RemoveDOFs = obj.LeftDOFsY;
            else
                error("Give a non-zero vector for the normals!")
            end
            % all others we keep
            obj.KeepDOFs = setdiff((1:obj.ParentElement.NumberDOF), obj.RemoveDOFs);
            obj.Name = "PeriodicNormal";

            Cleft = sparse(1:numel(obj.LeftDOFsX), obj.LeftDOFsX, nL(1), numel(obj.LeftDOFsX), parentElement.NumberDOF) + ...
                    sparse(1:numel(obj.LeftDOFsY), obj.LeftDOFsY, nL(2), numel(obj.LeftDOFsY), parentElement.NumberDOF);
            Cright = sparse(1:numel(obj.RightDOFsX), obj.RightDOFsX, nR(1), numel(obj.RightDOFsX), parentElement.NumberDOF) + ...
                     sparse(1:numel(obj.RightDOFsY), obj.RightDOFsY, nR(2), numel(obj.RightDOFsY), parentElement.NumberDOF);
            obj.C = Cleft + Cright;

            obj.CK = obj.C(:, obj.KeepDOFs);
            obj.CR = obj.C(:, obj.RemoveDOFs);
        end

        function [b, db]  = getBCvalues(obj, t)
            b = sparse(numel(obj.RemoveDOFs), 1);
            db = sparse(numel(obj.RemoveDOFs), 1);
        end
    end
end