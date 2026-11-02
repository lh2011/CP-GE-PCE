classdef IGA_BCdirichlet < IGA_BC & BC_Dirichlet
    properties (Access = public)
        Spaces, SpacesEval;
        Meshes, MeshesEval;
        const_f;
    end

    methods (Access = public)
        function obj = IGA_BCdirichlet(parentElement, boundaryNumbers, args)
            arguments
                parentElement IGA_Region
                boundaryNumbers 
                % for transient problems, derivative needs to be given: 
                % value = @(t) @(x,y) function of size x, depending of t, x, y
                % dvalue = @(t) @(x,y) function of size x, depending of t, x, y
                % for static problems:
                % value = @(x,y) function of size x
                % value = const
                args.value = 0;
                args.dvalue = 0;
            end
            boundaryPatchNumbers = parentElement.getBoundaryPatchNumbers(boundaryNumbers);
            dDOFs = parentElement.getBoundaryDOFs(boundaryPatchNumbers);
            
            obj@BC_Dirichlet(parentElement, dDOFs, value=args.value, dvalue=args.dvalue);
            obj.BoundaryPatchNumbers = boundaryPatchNumbers;
            obj.PlotColor = TUDa_getColor("2b");
            
            % Generate Spaces and Meshes
            for iBnd = 1:numel(boundaryPatchNumbers)
                patch = parentElement.PatchBoundaries(boundaryPatchNumbers(iBnd)).patches;
                side = parentElement.PatchBoundaries(boundaryPatchNumbers(iBnd)).faces;
                obj.MeshesEval{iBnd} = msh_eval_boundary_side(parentElement.Meshes.msh_patch{patch}, side);
                msh_side_int = msh_boundary_side_from_interior(parentElement.Meshes.msh_patch{patch}, side);
                obj.Spaces{iBnd} = parentElement.Spaces.sp_patch{patch}.constructor(msh_side_int);
                obj.SpacesEval{iBnd} = sp_precompute(obj.Spaces{iBnd}, msh_side_int, 'value', true, 'gradient', true);
            end

            % Generate IGA projection matrix:
            M = sparse(parentElement.NumberDOF, parentElement.NumberDOF);
            obj.const_f = sparse(parentElement.NumberDOF, 1);       % precalculate integral for spatially constant values
            for iBnd = 1:numel(boundaryPatchNumbers)
                patch = parentElement.PatchBoundaries(boundaryPatchNumbers(iBnd)).patches;
                patchDOFs = parentElement.Spaces.gnum{patch};
                M(patchDOFs, patchDOFs) = M(patchDOFs, patchDOFs) + ...
                    + parentElement.Length*int_u_v (obj.SpacesEval{iBnd}, obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, 1);
                obj.const_f(patchDOFs) = obj.const_f(patchDOFs) + parentElement.Length*int_f_v(obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd});
            end
            
            obj.CR = M(obj.RemoveDOFs, obj.RemoveDOFs);
            obj.CK = sparse(numel(obj.RemoveDOFs), numel(obj.KeepDOFs));
        end

        function [b, db] = getBCvalues(obj, t)
            if isa(obj.ValueDerivative, "function_handle")
               % Transient problem
               func = obj.Value(t);
               funcDer = obj.ValueDerivative(t);
               if isa(func, "function_handle")
                    % Spatially changing value:
                    b = sparse(obj.ParentElement.NumberDOF, 1);
                    for iBnd = 1:numel(obj.BoundaryPatchNumbers)
                        patch = obj.ParentElement.PatchBoundaries(obj.BoundaryPatchNumbers(iBnd)).patches;
                        patchDOFs = obj.ParentElement.Spaces.gnum{patch};
                        b(patchDOFs) = b(patchDOFs) + obj.ParentElement.Length*int_f_v(obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, func);
                    end
                    b = b(obj.RemoveDOFs);
                    if nargout == 2
                        db = sparse(obj.ParentElement.NumberDOF, 1);
                        for iBnd = 1:numel(obj.BoundaryPatchNumbers)
                            patch = obj.ParentElement.PatchBoundaries(obj.BoundaryPatchNumbers(iBnd)).patches;
                            patchDOFs = obj.ParentElement.Spaces.gnum{patch};
                            db(patchDOFs) = db(patchDOFs) + obj.ParentElement.Length*int_f_v(obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, funcDer);
                        end
                        db = db(obj.RemoveDOFs);
                    end
                else
                    % Spatially constant value:
                    b = func*obj.const_f(obj.RemoveDOFs);
                    if nargout == 2
                        db = funcDer*obj.const_f(obj.RemoveDOFs);
                    end
               end
            else
                % Static problem                
                if isa(obj.Value, "function_handle")
                    % Spatially changing value:
                    b = sparse(obj.ParentElement.NumberDOF, 1);
                    for iBnd = 1:numel(obj.BoundaryPatchNumbers)
                        patch = obj.ParentElement.PatchBoundaries(obj.BoundaryPatchNumbers(iBnd)).patches;
                        patchDOFs = obj.ParentElement.Spaces.gnum{patch};
                        b(patchDOFs) = b(patchDOFs) + obj.ParentElement.Length*int_f_v(obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, obj.Value);
                    end
                    b = b(obj.RemoveDOFs);
                    if nargout == 2
                        db = obj.ValueDerivative*obj.const_f(obj.RemoveDOFs);
                    end
                else
                    % Spatially constant value:
                    b = obj.Value*obj.const_f(obj.RemoveDOFs);
                    if nargout == 2
                        db = obj.ValueDerivative*obj.const_f(obj.RemoveDOFs);
                    end
                end
            end
        end
    end
end