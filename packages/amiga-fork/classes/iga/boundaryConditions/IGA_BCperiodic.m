classdef IGA_BCperiodic < IGA_BC & BC_Periodic

    methods (Access = public)
        function obj = IGA_BCperiodic(parentElement, sidesOne, sidesTwo)
            
            bndPatchesOne = parentElement.getBoundaryPatchNumbers(sidesOne);
            bndPatchesTwo = parentElement.getBoundaryPatchNumbers(sidesTwo);
            [sideOneDOFs, sideTwoDOFs] = parentElement.getBoundaryDOFsperiodic(bndPatchesOne, bndPatchesTwo);

            obj@BC_Periodic(parentElement, sideOneDOFs, sideTwoDOFs);

            obj.BoundaryPatchNumbers = union(bndPatchesOne, bndPatchesTwo);
            obj.PlotColor = TUDa_getColor("7b");
        end
    end
end