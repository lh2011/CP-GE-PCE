classdef BC_AntiPeriodic < BC
    properties (Access = public)
        SideOneDOFs, 
        SideTwoDOFs;
    end

    methods (Access = public)
        function obj = BC_AntiPeriodic(parentElement, sideOneDOFs, sideTwoDOFs)
            sideOneDOFs = reshape(sideOneDOFs, 1, []); sideTwoDOFs = reshape(sideTwoDOFs, 1, []);
            assert(isempty(intersect(sideOneDOFs, sideTwoDOFs)), "DOFs of both sides must not intersect!");
            % Reduce DOFs for possible multiplicity with other Boundary Conditions
            obj@BC(parentElement, union(sideOneDOFs, sideTwoDOFs, 'stable'));
            obj.SideOneDOFs = intersect(sideOneDOFs, obj.BoundaryDOFs, 'stable');
            obj.SideTwoDOFs = intersect(sideTwoDOFs, obj.BoundaryDOFs, 'stable');
            % only set SideOneDOFs as the DOFs which should be removed
            obj.RemoveDOFs = obj.SideOneDOFs;
            % all others are independent
            obj.KeepDOFs = setdiff((1:obj.ParentElement.NumberDOF), obj.RemoveDOFs);
            obj.Name = "AntiPeriodic";

            C = sparse(repmat(1:numel(obj.SideOneDOFs), 1, 2), [obj.SideOneDOFs, obj.SideTwoDOFs], ...
                ones(1, numel(obj.SideOneDOFs)+numel(obj.SideTwoDOFs)), ...
                numel(obj.RemoveDOFs), numel(obj.RemoveDOFs)+numel(obj.KeepDOFs));
            obj.CK = C(:, obj.KeepDOFs);
            obj.CR = C(:, obj.RemoveDOFs);
        end

        function [b, db]  = getBCvalues(obj, t)
            b = sparse(numel(obj.RemoveDOFs), 1);
            db = sparse(numel(obj.RemoveDOFs), 1);
        end
    end
end