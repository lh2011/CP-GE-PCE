% Right now only for spatially constant values!

classdef FEM_BCdirichlet < FEM_BC & BC_Dirichlet
    properties (Access = public)
        
    end

    methods (Access = public)
        function obj = FEM_BCdirichlet(parentElement, boundaryNumbers, args)
            arguments
                parentElement FEM_Region;
                boundaryNumbers 
                args.value = 0;     %  spatially constant or time dependent value
                args.dvalue = 0;    % spatially constant or time derivative of value
            end
            dDOFs = parentElement.getBoundaryDOFs(boundaryNumbers);
            
            obj@BC_Dirichlet(parentElement, dDOFs, value=args.value, dvalue=args.dvalue);
            obj.BoundaryNumbers = boundaryNumbers;
            obj.PlotColor = TUDa_getColor("2b");
        end
    end
end