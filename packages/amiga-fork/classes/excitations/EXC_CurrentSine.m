classdef EXC_CurrentSine < EXC_Current
    properties (Access = public)
        IApp;
        Omega;
        UVWshift;
        Phase0;
    end

    methods
        function obj = EXC_CurrentSine(material, phasesStruct, iPeak, omega, uvwshift, phase0, args)
            arguments
                material 
                phasesStruct 
                iPeak 
                omega 
                uvwshift 
                phase0
                args.NumberWindings = 1;
                args.PlotColor = [0,0,0];
            end
            currentFunction = @(t) iPeak*sin(omega*t + uvwshift + phase0);
            obj@EXC_Current(material, phasesStruct, currentFunction, "NumberWindings", args.NumberWindings, "PlotColor", args.PlotColor)
            obj.IApp = iPeak;
            obj.Omega = omega;
            obj.UVWshift = uvwshift;
            obj.Phase0 = phase0;
            % Update current function with object values
            obj.ExcitationFunction = @(t) obj.IApp*sin(obj.Omega*t + obj.UVWshift + obj.Phase0);
        end
    end
end