function materials = parseFEMMthermal(optionalName)
    % Open the file
    filename = 'heatlib.dat';
    fid = fopen(filename, 'r');
    if fid == -1
        error('Failed to open file.');
    end

    materials = [];
    currentMaterial = struct();
    insideBlock = false;
    readTK = false;
    tkCounter = 0;
    tkData = [];

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        
        % Start of a new block
        if contains(line, '<BeginBlock>')
            insideBlock = true;
            currentMaterial = struct();
            tkData = [];
            continue;
        end

        % End of a block
        if contains(line, '<EndBlock>') && insideBlock
            if isfield(currentMaterial, 'BlockName')
                if ~isempty(tkData)
                    currentMaterial.TK = tkData;
                else
                    currentMaterial.TK = [];
                    currentMaterial.TKPoints = 0;
                end
                if isempty(materials)
                    materials = currentMaterial;
                else
                    materials(end+1) = currentMaterial; %#ok<AGROW>
                end
            end
            insideBlock = false;
            readTK = false;
            continue;
        end

        if insideBlock && ~readTK
            % Parse scalar values
            tokens = regexp(line, '<(?<key>\w+)> = (?<value>.+)', 'names');
            if ~isempty(tokens)
                key = tokens.key;
                valStr = strrep(tokens.value, '"', ''); % remove quotes
                val = str2double(valStr);
                if isnan(val)
                    currentMaterial.(key) = valStr; % keep as string
                else
                    currentMaterial.(key) = val;    % convert to number
                end
                % check if only a certain material is wanted
                if nargin > 0 && isfield(currentMaterial, "BlockName")
                    if ~strcmp(optionalName, currentMaterial.BlockName)
                        insideBlock = false;
                        continue
                    end
                end

                if strcmp(key, 'TKPoints')
                    tkCounter = str2double(valStr);
                    if tkCounter > 0
                        readTK = true;
                        tkData = zeros(tkCounter, 2);
                        tkLineIndex = 1;
                    end
                end
            end

        elseif readTK
            % Read TK data lines
            if tkLineIndex <= tkCounter
                tkValues = sscanf(line, '%f');
                if length(tkValues) == 2
                    tkData(tkLineIndex, :) = tkValues';
                end
                tkLineIndex = tkLineIndex + 1;
            end
            if tkLineIndex > tkCounter
                readTK = false;
            end
        end
    end

    fclose(fid);
end

