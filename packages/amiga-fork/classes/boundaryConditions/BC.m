classdef BC < handle
    properties (Access = public)
        Name;
        ParentElement;      % DAE_Element to apply the BC
        BoundaryDOFs;       % Alls DOFs of the boundary
        RemoveDOFs;         % DOFs which will be removed from the equation system
        KeepDOFs;           % DOFs which remain in the equation system
        C;                  % General boundary condition for C*y = b
        CK;                 % C(:, KeepDOFs)
        CR                  % C(:, RemoveDOFs)
    end

    methods (Abstract)
        [b, db]  = getBCvalues(obj, t)
    end

    methods (Access = public)
        function obj = BC(parentElement, boundaryDOFs)
            obj.ParentElement = parentElement;
            obj.BoundaryDOFs = reshape(boundaryDOFs, 1, []);
            % Check all other boundary conditions for consistency
            for iBC = 1:numel(obj.ParentElement.BoundaryConditions)
                checkBC = obj.ParentElement.BoundaryConditions(iBC).BoundaryCondition;
                % occupiedDOFs = intersect(obj.BoundaryDOFs, checkBC.RemoveDOFs, 'stable'); % TBD: return CK AND CR for this procedure, otherwise not correct G/H
                occupiedDOFs = intersect(obj.BoundaryDOFs, checkBC.BoundaryDOFs, 'stable');
                if ~isempty(occupiedDOFs)
                    obj.BoundaryDOFs = setdiff(obj.BoundaryDOFs, occupiedDOFs, 'stable');
                    disp(['Note: DOF(s) ' num2str(reshape(occupiedDOFs, 1, [])) ' were ignored in setting the boundary condition since they are alreay in use.'])
                end
            end
            parentElement.addBoundaryCondition(obj);
        end

        function [dofsD, dofsI]  = getBCindices(obj)
            dofsD = obj.RemoveDOFs;
            dofsI = obj.KeepDOFs;
        end

        function [G, H] = getBCmatrices(obj)
            H = inv(obj.CR);
            G = -H*obj.CK;
        end
    end
end