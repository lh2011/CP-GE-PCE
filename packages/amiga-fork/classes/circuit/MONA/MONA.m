classdef  MONA < EQS_System
    properties
        Nodes
    end
    
    methods 
        function obj = MONA()
            obj.NumberDOF = 0;
            obj.Nodes = MONA_Nodes();
            obj.Elements = struct('Element', {obj.Nodes}, 'Indices', {[]}, 'IndicesExt', {[]}, 'IndicesInt', {[]});
        end

        % Extend addElement function of EQS_System for the nodes
        function addElement(obj, element, nodes)
            arguments
                obj 
                element {mustBeA(element, 'MONA_Element')};
                nodes 
            end
            idxNodes = obj.addNodes(nodes);
            element.addExternalElement(obj.Nodes, idxNodes)

            addElement@EQS_System(obj, element);
            % Set element name
            if isempty(element.Name)
                numElements = sum(cellfun(@(c) isa(c, class(element)), {obj.Elements.Element}));
                element.Name = [num2str(class(element)) '_' num2str(numElements)];
            end
        end

        function idxNodes = addNodes(obj, nodes)
            [idxNodes, idxNewNodes] = obj.Nodes.addNodes(nodes);
            obj.Elements(1).IndicesInt = [obj.Elements(1).IndicesInt, obj.NumberDOF + (1:numel(idxNewNodes))];
            obj.Elements(1).Indices = obj.Elements(1).IndicesInt;
            obj.NumberDOF = obj.NumberDOF + numel(idxNewNodes);
        end
    end
end


