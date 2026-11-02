classdef IGA_BCperiodicNormal < IGA_BC & BC_PeriodicNormal

    methods (Access = public)
        function obj = IGA_BCperiodicNormal(parentElement, leftSides, rightSides, nLeft, nRight)

            boundaryPatchNumbersLeft = parentElement.getBoundaryPatchNumbers(leftSides);
            boundaryPatchNumbersRight = parentElement.getBoundaryPatchNumbers(rightSides);

            [leftDOFsX, leftDOFsY, rightDOFsX, rightDOFsY]  = parentElement.getBoundaryDOFsperiodic(boundaryPatchNumbersLeft, boundaryPatchNumbersRight);

            obj@BC_PeriodicNormal(parentElement, leftDOFsX, leftDOFsY, rightDOFsX, rightDOFsY, nLeft, nRight);

            obj.BoundaryPatchNumbers = union(boundaryPatchNumbersLeft, boundaryPatchNumbersRight);
            obj.PlotColor = TUDa_getColor("11c");
        end
    end
end

