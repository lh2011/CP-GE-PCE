classdef FEM_BCfloating < FEM_BC & BC_Floating

    methods (Access = public)
        function obj = FEM_BCfloating(parentElement, boundaryNumbers)
            
            sideDOFs = parentElement.getBoundaryDOFs(boundaryNumbers);

            obj@BC_Floating(parentElement, sideDOFs);

            obj.BoundaryNumbers = boundaryNumbers;
            obj.PlotColor = TUDa_getColor("10a");
        end
    end
end