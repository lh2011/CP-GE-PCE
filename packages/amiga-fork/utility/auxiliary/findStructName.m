function index = findStructName(myStruct, name)
    index = [];
    for iStruct = 1:numel(myStruct)
        if strcmp(myStruct(iStruct).("Name"), name)
            index = iStruct;
            break;
        end
    end
end