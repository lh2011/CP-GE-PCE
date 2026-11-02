classdef  MLA < EQS_System
    properties
        DirectedGraph=digraph();
        UndirectedGraph;
        Loops;
    end
    
    methods 
        function obj = MLA()
            obj.NumberDOF = 0;
            obj.Loops = MLA_Loops();

            obj.Elements = struct('Element', {obj.Loops}, 'Indices', {[]}, 'IndicesExt', {[]}, 'IndicesInt', {[]});
        end

        % Extend addElement function of EQS_System for the nodes
        function addElement(obj, element, nodes, updateLoops)
            arguments
                obj 
                element 
                nodes 
                updateLoops (1,1) logical = true;
            end

            % Set element name
            if isempty(element.Name)
                numElements = sum(cellfun(@(c) isa(c, class(element)), {obj.Elements.Element}));
                element.Name = [num2str(class(element)) '_' num2str(numElements+1)];
            end
            if ~strcmp(element.Name(1:3), 'MLA')
                warning('MLA:Element', 'You might have added a non MLA element to the MLA!');
            end

            elementIndex = numel(obj.Elements) + 1;
            elementTable = table([string(nodes(1)), string(nodes(2))], string(element.Name), elementIndex, 'VariableNames', {'EndNodes', 'ElementName', 'ElementIndex'});
            obj.DirectedGraph = obj.DirectedGraph.addedge(elementTable);

            % add new element
            element.addExternalElement(obj.Loops, [])
            addElement@EQS_System(obj, element);

            if updateLoops
                obj.updateLoops()
            end
        end

        function updateLoops(obj)
            elementTable = table(obj.DirectedGraph.Edges.EndNodes, obj.DirectedGraph.Edges.ElementName, obj.DirectedGraph.Edges.ElementIndex, 'VariableNames', {'EndNodes', 'ElementName', 'ElementIndex'});
            obj.UndirectedGraph = graph(elementTable);
            % edge sorting directed graph to undirected graph
            [~, u2d] = ismember(obj.UndirectedGraph.Edges.ElementIndex, obj.DirectedGraph.Edges.ElementIndex);
            [~, u2el] = ismember(obj.UndirectedGraph.Edges.ElementIndex, sort(obj.UndirectedGraph.Edges.ElementIndex));
            % nodeReversedUD = vertcat(obj.DirectedGraph.Edges.EndNodes{u2d,1}) ~= vertcat(obj.UndirectedGraph.Edges.EndNodes{:, 1});
            % sanityCheck = vertcat(obj.DirectedGraph.Edges.EndNodes{u2d,1}) ~= vertcat(obj.UndirectedGraph.Edges.EndNodes{:, 2});
            nodeReversedUD = cellfun(@strcmp, obj.DirectedGraph.Edges.EndNodes(u2d,1), obj.UndirectedGraph.Edges.EndNodes(:, 1));
            sanityCheck = cellfun(@strcmp, obj.DirectedGraph.Edges.EndNodes(u2d,1), obj.UndirectedGraph.Edges.EndNodes(:, 2));
            assert(all(nodeReversedUD~=sanityCheck));
            % check that nodes are not reordered
            % assert(all(vertcat(obj.DirectedGraph.Nodes.Name{:})==vertcat(obj.UndirectedGraph.Nodes.Name{:})));
            assert(all(cellfun(@strcmp, obj.DirectedGraph.Nodes.Name, obj.UndirectedGraph.Nodes.Name)));

            [cycles, edgecycles] = cyclebasis(obj.UndirectedGraph);
            % Update Loops
            obj.Loops.NumberDOF = numel(cycles);
            obj.Elements(1).IndicesInt = 1:obj.Loops.NumberDOF;
            obj.Elements(1).Indices = 1:obj.Loops.NumberDOF;
            % Reset Element Vectors and external dependencies
            for element = elementTable.ElementIndex'
                obj.Elements(element).Element.B = [];
                obj.Elements(element).Element.NumberDOFext = 0;
                obj.Elements(element).Indices = [];
                obj.Elements(element).IndicesExt = [];
                obj.Elements(element).IndicesInt = [];
            end
            % Update Elements in Loops
            for iCycle = 1:numel(cycles)
                cycle = cycles{iCycle};
                edges = edgecycles{iCycle};
                if isscalar(edges)
                    error("You have an element connected to itself");
                end
                edgeSign = zeros(numel(edges), 1);
                % first edge defines direction
                edgeSign(1) = 1;
                % sign for each edge in the cycle
                for iEdge = 2:numel(edges)
                    prevEdge = edges(iEdge-1);
                    edge = edges(iEdge);
                    % find correct sign of cycle
                    % edgeReversed = (obj.UndirectedGraph.Edges.EndNodes{edge,1} == obj.UndirectedGraph.Edges.EndNodes{prevEdge,1}) ||  (obj.UndirectedGraph.Edges.EndNodes{edge,2} == obj.UndirectedGraph.Edges.EndNodes{prevEdge,2});
                    edgeReversed = strcmp(obj.UndirectedGraph.Edges.EndNodes{edge,1}, obj.UndirectedGraph.Edges.EndNodes{prevEdge,1}) ||  strcmp(obj.UndirectedGraph.Edges.EndNodes{edge,2}, obj.UndirectedGraph.Edges.EndNodes{prevEdge,2});
                    if edgeReversed
                        edgeSign(iEdge) = -edgeSign(iEdge-1);
                    else
                        edgeSign(iEdge) = edgeSign(iEdge-1);
                    end
                end
                % check if directed graph is differently ordered than undirected graph
                d2uSign = -2*nodeReversedUD(edges) +1;
                edgeSign = edgeSign .* d2uSign;
                
                % update indices for elements
                for iEdge = 1:numel(edges)
                    idxEl = u2el(edges(iEdge)) + 1;
                    % increase number of DOF because element is in this loop
                    obj.Elements(idxEl).Element.NumberDOFext = obj.Elements(idxEl).Element.NumberDOFext + 1;
                    % update external dependencies of elements to loop
                    obj.Elements(idxEl).IndicesExt = [obj.Elements(idxEl).IndicesExt, iCycle];
                    % update B vector
                    obj.Elements(idxEl).Element.B = [obj.Elements(idxEl).Element.B; edgeSign(iEdge)];
                end
            end
            % update internal indices for all elements
            obj.NumberDOF = 0;
            for iEl = 1:numel(obj.Elements)
                obj.Elements(iEl).IndicesInt = obj.NumberDOF+(1:obj.Elements(iEl).Element.NumberDOF);
                obj.NumberDOF = obj.NumberDOF + obj.Elements(iEl).Element.NumberDOF;
                obj.Elements(iEl).Indices = [obj.Elements(iEl).IndicesExt, obj.Elements(iEl).IndicesInt];
            end
        end

        function plotGraph(obj)
            figure
            obj.DirectedGraph.plot('EdgeLabel', obj.DirectedGraph.Edges.ElementName, 'interpreter', 'none');
        end

        function plotLoops(obj)
            figure
            [cycles,edgecycles] = cyclebasis(obj.UndirectedGraph);
            tiledlayout flow
            for k = 1:length(cycles)
                nexttile
                highlight(plot(obj.UndirectedGraph),cycles{k},'Edges',edgecycles{k},'EdgeColor','r','NodeColor','r')
                title("Cycle " + k)
            end
        end
    end
end


