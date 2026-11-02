classdef FEM_BCperiodic < FEM_BC & BC_Periodic

    methods (Access = public)
        function obj = FEM_BCperiodic(parentElement, sidesOne, sidesTwo)
            
            [leftDOFs, rightDOFs] = parentElement.getBoundaryDOFsperiodic(sidesOne, sidesTwo);

            obj@BC_Periodic(parentElement, leftDOFs, rightDOFs);

            obj.BoundaryNumbers = union(sidesOne, sidesTwo);
            obj.PlotColor = TUDa_getColor("7b");
        end
    end
end