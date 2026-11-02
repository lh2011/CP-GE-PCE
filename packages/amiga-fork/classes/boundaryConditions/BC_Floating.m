classdef BC_Floating < BC
    properties (Access = public)
        
    end

    methods (Access = public)
        function obj = BC_Floating(parentElement, boundaryDOFs)
            obj@BC(parentElement, boundaryDOFs);
            obj.RemoveDOFs = obj.BoundaryDOFs(2:end);
            obj.KeepDOFs = setdiff(1:parentElement.NumberDOF, obj.RemoveDOFs);

            obj.Name = "Floating";

            C = sparse(numel(obj.BoundaryDOFs)-1, parentElement.NumberDOF);
            C(:, obj.BoundaryDOFs(1:end-1)) = diag(ones(numel(obj.BoundaryDOFs)-1, 1));
            C(:, obj.BoundaryDOFs(2:end)) = C(:, obj.BoundaryDOFs(2:end)) - diag(ones(numel(obj.BoundaryDOFs)-1, 1));

            obj.CK = C(:, obj.KeepDOFs);
            obj.CR = C(:, obj.RemoveDOFs);
        end

        function [b, db]  = getBCvalues(obj, t)
            b = zeros(numel(obj.RemoveDOFs),1);
            db = zeros(numel(obj.RemoveDOFs),1);
        end
    end
end