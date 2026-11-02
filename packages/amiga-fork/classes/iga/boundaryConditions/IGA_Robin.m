% Only constant values allowed for now in space and time

classdef IGA_Robin < EQS_Element
    properties (Access = public)
        BoundaryNumbers;
        PlotColor;
        MeshesEval;
        SpacesEval;
        Spaces;
        ParentElement;
        RefValue, Flux, Gamma;
    end

    methods (Access = public)
        function obj = IGA_Robin(parentElement, boundaryNumbers, refValue, flux, gamma)
            if ~exist("refValue", "var") || isempty(refValue)
                refValue = 0;
                warning('IGA_Robin:NoRefValue', 'No refrence value given in Robin BC. This corresponds to a Neumann Boundary Condition.\n')
            end
            obj.RefValue = refValue;
            if ~exist("flux", "var") || isempty(flux)
                flux = 0;
                warning('IGA_Robin:NoFlux', 'No flux given in Robin BC. This corresponds to a Robin Boundary Condition without flux.\n')
            end
            obj.Flux = flux;
            if ~exist("gamma", "var") || isempty(gamma)
                gamma = 0;
                warning('IGA_Robin:NoGamma', 'No gamma given in Robin BC. This corresponds to a Neumann Boundary Condition.\n')
            end
            obj.Gamma = gamma;

            obj.NumberDOF = 0;
            obj.ParentElement = parentElement;
            obj.addExternalElement(parentElement);
            
            obj.BoundaryNumbers = boundaryNumbers;
            obj.PlotColor = TUDa_getColor("10b");

            % Generate Spaces and Meshes
            for iBnd = 1:numel(obj.BoundaryNumbers)
                patch = obj.ParentElement.PatchBoundaries(obj.BoundaryNumbers(iBnd)).patches;
                side = obj.ParentElement.PatchBoundaries(obj.BoundaryNumbers(iBnd)).faces;
                obj.MeshesEval{iBnd} = msh_eval_boundary_side(obj.ParentElement.Meshes.msh_patch{patch}, side);
                msh_side_int = msh_boundary_side_from_interior(obj.ParentElement.Meshes.msh_patch{patch}, side);
                obj.Spaces{iBnd} = obj.ParentElement.Spaces.sp_patch{patch}.constructor(msh_side_int);
                obj.SpacesEval{iBnd} = sp_precompute(obj.Spaces{iBnd}, msh_side_int, 'value', true, 'gradient', false);
            end
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
            Mlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOFext, obj.NumberDOFext);
            if nargout == 2
                Klin = sparse(obj.NumberDOFext, obj.NumberDOFext);
                for iBnd = 1:numel(obj.BoundaryNumbers)
                    patch = obj.ParentElement.PatchBoundaries(obj.BoundaryNumbers(iBnd)).patches;
                    patchDOFs = obj.ParentElement.Spaces.gnum{patch};
                    Klin(patchDOFs, patchDOFs) = Klin(patchDOFs, patchDOFs) + ... 
                        - obj.ParentElement.Length*int_u_v(obj.SpacesEval{iBnd}, obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, obj.Gamma);
                end
            end
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext, obj.NumberDOFext);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOFext, 1);
            if nargout == 2
                Flin = sparse(obj.NumberDOFext, 1);
                for iBnd = 1:numel(obj.BoundaryNumbers)
                    patch = obj.ParentElement.PatchBoundaries(obj.BoundaryNumbers(iBnd)).patches;
                    patchDOFs = obj.ParentElement.Spaces.gnum{patch};
                    Flin(patchDOFs) = Flin(patchDOFs) + ...
                        obj.ParentElement.Length*int_f_v(obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, obj.Gamma*obj.RefValue + obj.Flux);
                end
            end
        end

        function postprocess(obj, t, ys, ysprime)
            % TBD energy flow: Int gamma * (Tref-T) + h dG
        end
    end
end