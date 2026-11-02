function materials = parseFEMMmagnetic(optionalName)
    % Open the file
    filename = 'matlib.dat';
    fid = fopen(filename, 'r');
    if fid == -1
        error('Failed to open file.');
    end

    materials = [];
    currentMaterial = struct();
    insideBlock = false;
    readBH = false;
    bhCounter = 0;
    bhData = [];

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        
        % Start of a new block
        if contains(line, '<BeginBlock>')
            insideBlock = true;
            currentMaterial = struct();
            bhData = [];
            continue;
        end

        % End of a block
        if contains(line, '<EndBlock>') && insideBlock
            if isfield(currentMaterial, 'BlockName')
                if ~isempty(bhData)
                    currentMaterial.BH = bhData;
                else
                    currentMaterial.BH = [];
                end
                if isempty(materials)
                    materials = currentMaterial;
                else
                    materials(end+1) = currentMaterial; %#ok<AGROW>
                end
            end
            insideBlock = false;
            readBH = false;
            continue;
        end

        if insideBlock && ~readBH
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

                if strcmp(key, 'BHPoints')
                    bhCounter = str2double(valStr);
                    if bhCounter > 0
                        readBH = true;
                        bhData = zeros(bhCounter, 2);
                        bhLineIndex = 1;
                    end
                end
            end

        elseif readBH
            % Read BH data lines
            if bhLineIndex <= bhCounter
                bhValues = sscanf(line, '%f');
                if length(bhValues) == 2
                    bhData(bhLineIndex, :) = bhValues';
                end
                bhLineIndex = bhLineIndex + 1;
            end
            if bhLineIndex > bhCounter
                readBH = false;
            end
        end
    end

    fclose(fid);
end

