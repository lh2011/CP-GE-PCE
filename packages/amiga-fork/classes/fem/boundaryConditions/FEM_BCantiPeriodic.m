classdef FEM_BCantiPeriodic < FEM_BC & BC_AntiPeriodic

    methods (Access = public)
        function obj = FEM_BCantiPeriodic(parentElement, sidesOne, sidesTwo)

            [leftDOFs, rightDOFs] = parentElement.getBoundaryDOFsperiodic(sidesOne, sidesTwo);

            obj@BC_AntiPeriodic(parentElement, leftDOFs, rightDOFs);

            obj.BoundaryNumbers = union(sidesOne, sidesTwo);
            obj.PlotColor = TUDa_getColor("5a");
        end
    end
end