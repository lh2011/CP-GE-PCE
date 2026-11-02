classdef IGA_BCfloating < IGA_BC & BC_Floating

    methods (Access = public)
        function obj = IGA_BCfloating(parentElement, boundaryNumbers)

            boundaryPatchNumbers = parentElement.getBoundaryPatchNumbers(boundaryNumbers);
            sideDOFs = parentElement.getBoundaryDOFs(boundaryPatchNumbers);

            obj@BC_Floating(parentElement, sideDOFs);

            obj.BoundaryPatchNumbers = boundaryPatchNumbers;
            obj.PlotColor = TUDa_getColor("10a");
        end
    end
end