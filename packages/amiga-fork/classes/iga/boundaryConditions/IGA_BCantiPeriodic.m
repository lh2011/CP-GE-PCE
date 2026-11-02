classdef IGA_BCantiPeriodic < IGA_BC & BC_AntiPeriodic

    methods (Access = public)
        function obj = IGA_BCantiPeriodic(parentElement, sidesOne, sidesTwo)

            bndPatchesOne = parentElement.getBoundaryPatchNumbers(sidesOne);
            bndPatchesTwo = parentElement.getBoundaryPatchNumbers(sidesTwo);
            [sideOneDOFs, sideTwoDOFs] = parentElement.getBoundaryDOFsperiodic(bndPatchesOne, bndPatchesTwo);

            obj@BC_AntiPeriodic(parentElement, sideOneDOFs, sideTwoDOFs);

            obj.BoundaryPatchNumbers = union(bndPatchesOne, bndPatchesTwo);
            obj.PlotColor = TUDa_getColor("5a");
        end
    end
end