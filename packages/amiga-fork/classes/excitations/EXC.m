classdef EXC < handle
    properties (Access = public)
        Name;
        ParentRegion;
        ParentMaterial;
        RegionNumbers;
        LocalIndices;
        Space, Mesh;
        ExcitationFunction;
        IsTimeDependent = false;
        IsSpaceDependent = false;
        RHS;                            % generic vector for constant excitations
    end


    methods (Access = public)
        function obj = EXC(material, regionNumbers)
            arguments
                material MAT
                regionNumbers = material.RegionNumbers;
            end
            assert(all(ismember(regionNumbers, material.RegionNumbers)), "You are setting the excitation for a region, which is not defined the given material!");
            obj.ParentMaterial = material;
            obj.ParentRegion = material.Region;
            obj.RegionNumbers = regionNumbers;
            [obj.Mesh, obj.Space, obj.LocalIndices] = material.extractMeshSpace(regionNumbers);
            material.addExcitation(obj);
        end
        
        % possible inputs:
        % valueOrFunction = value
        % valueOrFunction = @(t) value(t)
        % valueOrFunction = @(x,y) value(x,y)
        % valueOrFunction = @(t) @(x,y) value(t,x,y)
        function setValue(obj, valueOrFunction)
            % check input
            if isa(valueOrFunction, "function_handle")
                if nargin(valueOrFunction) == 1
                    obj.IsTimeDependent = true;
                    spaceFunction = valueOrFunction(0);
                    if isa(spaceFunction, "function_handle")
                        assert(nargin(spaceFunction) == 2, "function structure must be @(t) @(x,y) ...");
                        obj.IsSpaceDependent = true;
                    else
                        obj.IsSpaceDependent = false;
                    end
                elseif nargin(valueOrFunction) == 2
                    obj.IsTimeDependent = false;
                    obj.IsSpaceDependent = true;
                end
            else
                obj.IsTimeDependent = false;
                obj.IsSpaceDependent = false;
            end

            obj.ExcitationFunction = valueOrFunction;
            obj.RHS = obj.Mesh.length*int_f_v(obj.Space, obj.Mesh, 1);
        end

        function [Fnonlin, Flin] = getEXCvalues(obj, t)
            if obj.IsSpaceDependent && obj.IsTimeDependent
                valueOfT = obj.ExcitationFunction(t);
                Fnonlin = obj.Mesh.length*int_f_v(obj.Space, obj.Mesh, valueOfT);
                if nargout == 2
                    Flin = zeros(obj.Space.ndof, 1);
                end
            elseif obj.IsSpaceDependent
                Fnonlin = zeros(obj.Space.ndof, 1);
                if nargout == 2
                    Flin = obj.Mesh.length*int_f_v(obj.Space, obj.Mesh, obj.ExcitationFunction);
                end
            elseif obj.IsTimeDependent
                value = obj.ExcitationFunction(t);
                Fnonlin = value*obj.RHS;
                if nargout == 2
                    Flin = zeros(obj.Space.ndof, 1);
                end
            else
                Fnonlin = zeros(obj.Space.ndof, 1);
                if nargout == 2
                    Flin = obj.ExcitationFunction*obj.RHS;
                end
            end
        end
    end
end