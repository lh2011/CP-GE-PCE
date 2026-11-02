classdef MAT < handle
    properties
        Rho = 0;                   % Material Density
        CostPerKG = 0;             % Material Cost
        PlotColor = [0.9, 0.9, 0.9];      
        PlotAlpha = 1;

        Region, RegionNumbers;
        Space, Mesh;
        Excitations = struct('Excitation', {}, 'Indices', {});
    end

    methods (Access = public)
        function rho = getRho(obj)
            rho = obj.Rho;
        end
        
        function rho = setRho(obj, rho)
            obj.Rho = rho;
        end

        function cost = getCostperKG(obj)
            cost = obj.CostPerKG;
        end

        function setCostperKG(obj, cost)
            obj.CostPerKG = cost;
        end

        function col = getPlotColor(obj)
            col = obj.PlotColor;
        end

        function alpha = getPlotAlpha(obj)
            alpha = obj.PlotAlpha;
        end

        function area = getArea(obj)
            area = int_Omega(obj.Mesh);
        end

        function volume = getVolume(obj)
            volume = obj.Mesh.length*obj.getArea();
        end

        function weight = getWeight(obj)
            weight = obj.Rho*obj.getVolume();
        end

        function cost = getCost(obj)
            cost = obj.CostPerKG*obj.getWeight();
        end

        function inertia = getInertia(obj)
            inertia = obj.Rho*obj.Mesh.length*sum(obj.Mesh.jacdet.*obj.Mesh.quad_weights .* reshape((obj.Mesh.geo_map(1, :, :).^2 + obj.Mesh.geo_map(2, :, :).^2) , obj.Mesh.nqn, obj.Mesh.nel), 'all');
        end

        function addExcitation(obj, excitation)
            pos = numel(obj.Excitations)+1;
            excitation.Name = inputname(2);
            obj.Excitations(pos).Excitation = excitation;
        end

        function resetExcitations(obj)
            obj.Excitations = struct('Excitation', {}, 'Indices', {});
            obj.Region.Klin = [];
            obj.Region.Flin = [];
        end

        % function for extracting a subspace of the mesh/space
        function [submesh, subspace, indices] = extractMeshSpace(obj, regionNumbers)
            assert(all(ismember(regionNumbers, obj.RegionNumbers)), "Can only extract subspaces of regions that belong to the material!");
            subRegions = obj.Region.getPatchNumbers(regionNumbers);
            submesh = MESH_2D(obj.Region.MeshesEval(subRegions), "length", obj.Region.Length);
            if nargout > 1
                subspace = SPACE_2D(obj.Region.SpacesEval(subRegions), obj.Region.Spaces.gnum(subRegions), obj.Region.Spaces.ncomp);
            end
            if nargout > 2
                [~, ~, indices] = intersect(subspace.gnum, obj.Space.gnum);
            end
        end
    end
end