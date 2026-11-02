% 

function params = TIND_rotor(params)
    arguments
        % Standard parameters   
        params.Symmetry = 4;
        params.BARS = 40;              % Number of cooling pipes
        params.MeshSize = 0;            % Automesh
        params.MeshSizeAG = 1e-3;        % mesh sides at air gap
    
        params.HB = 22.3e-3;            % Height Bar
        params.HSO = 4.0e-3;            % Height Slot Opening

        params.RAG = 186.75e-3;         % Radius Air Gap
        params.RPosPipe1 = 95e-3;       % Radius to center of pipes
        params.RPosPipe2 = 112.5e-3;    % Radius to center of pipes
        params.RPipe = 7.5e-3;          % Radius of cooling pipes
        params.RRI = 55e-3;             % Rotor Radius Inner
        params.RRO = 186e-3;            % Rotor Radius Outer
        params.WB = 12.3e-3;            % Width Bar
        params.WSO = 2.2e-3;            % Width Slot Opening
    end
    
    % calculated parameters
    a1 = asin((params.WSO/2)/(params.RRO));
    a2 = pi/params.BARS;
    a3 = a2-a1;
    barAngle = 2*pi/params.BARS;
    rotorAngle = barAngle*params.BARS/params.Symmetry;

    point{1} = params.RAG*[1, 0];
    point{2} = params.RRO*[1, 0]; 
    point{3} = params.RRI*[1, 0];
    point{4} = params.RAG*[cos(rotorAngle), sin(rotorAngle)];
    point{5} = params.RRO*[cos(rotorAngle), sin(rotorAngle)];
    point{6} = params.RRI*[cos(rotorAngle), sin(rotorAngle)];

    point{7} =  [cos(a1)*params.RRO, -params.WSO/2];
    point{8} =  [point{7}(1)-params.HSO, point{7}(2)];
    point{9} =  [point{8}(1), -params.WB/2];
    point{10} =  point{9} - [params.HB, 0];
    point{11} =  point{10} + [0, params.WB];
    point{12} =  point{9} + [0, params.WB];
    point{13} =  point{8} + [0, params.WSO];
    point{14} =  [point{7}(1), - point{7}(2)];
    point{15} = params.RRO*[cos(barAngle), sin(barAngle)];

    % correct unrotated points
    for iPoint = 7:15
        point{iPoint} = ([cos(barAngle/2), -sin(barAngle/2); sin(barAngle/2), cos(barAngle/2)]*point{iPoint}')';
    end
    nPoints = numel(point);
    % copy and rotate points for all BARS
    for iSlot = 1:(params.BARS /params.Symmetry)
        for iPoint = 1:nPoints
            rotAngle = barAngle*iSlot;
            point{nPoints+(iSlot-1)*nPoints+iPoint} =  ([cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)]*point{iPoint}')';
        end
    end

    % figure
    % hold on
    % for i = 1:numel(point)
    %     scatter(point{i}(1), point{i}(2), "black", "filled");
    %     text(point{i}(1)+0.0001, point{i}(2), string(i))
    % end

    openfemm();
    newdocument(0);
    mi_probdef(0, 'meters', 'planar' ,1e-6, 0.1 ,30, 0);

    % GENERATE GEOMETRY
    % outer line
    mi_drawarc(point{1}, point{4}, rad2deg(rotorAngle), 360*params.MeshSizeAG/params.RAG);
    mi_drawarc(point{3}, point{6}, rad2deg(rotorAngle), 1);

    mi_drawline(point{1}, point{2});
    mi_drawline(point{2}, point{3});
    mi_drawline(point{4}, point{5});
    mi_drawline(point{5}, point{6});

    mi_getmaterial('10 AWG');

    % create BARS
    for iSlot = 1:(params.BARS/params.Symmetry)
        mi_drawarc(point{2+(iSlot-1)*nPoints}, point{7+(iSlot-1)*nPoints}, rad2deg(a3), 360*params.MeshSizeAG/params.RRO);
        mi_drawline(point{7+(iSlot-1)*nPoints}, point{8+(iSlot-1)*nPoints});
        mi_drawline(point{8+(iSlot-1)*nPoints}, point{9+(iSlot-1)*nPoints});
        mi_drawline(point{9+(iSlot-1)*nPoints}, point{10+(iSlot-1)*nPoints});
        mi_drawline(point{10+(iSlot-1)*nPoints}, point{11+(iSlot-1)*nPoints});
        mi_drawline(point{11+(iSlot-1)*nPoints}, point{12+(iSlot-1)*nPoints});
        mi_drawline(point{12+(iSlot-1)*nPoints}, point{13+(iSlot-1)*nPoints});
        mi_drawline(point{13+(iSlot-1)*nPoints}, point{14+(iSlot-1)*nPoints});
        mi_drawline(point{8+(iSlot-1)*nPoints}, point{13+(iSlot-1)*nPoints});
        mi_drawarc(point{14+(iSlot-1)*nPoints}, point{17+(iSlot-1)*nPoints}, rad2deg(a3), 360*params.MeshSizeAG/params.RRO);
        % set copper
        xm = (point{10+(iSlot-1)*nPoints} + point{12+(iSlot-1)*nPoints})/2;
        mi_addblocklabel(xm);
        mi_selectlabel(xm);
        mi_setblockprop('10 AWG',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
        mi_clearselected();
    end

    % create cooling pipes, same number as bars
    nBars = floor(params.BARS/params.Symmetry) + 1;
    for iPipe = 1:nBars
        rotAngle = (iPipe-1) * 2*pi/params.BARS;
        if mod(iPipe, 2) == 0
            p1 = ([cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)]*[params.RPosPipe1 + params.RPipe; 0])';
            p2 = ([cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)]*[params.RPosPipe1 - params.RPipe; 0])';
        else
            p1 = ([cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)]*[params.RPosPipe2 + params.RPipe; 0])';
            p2 = ([cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)]*[params.RPosPipe2 - params.RPipe; 0])';
        end
        if iPipe == 1
            mi_drawarc(p1, p2, 180, 20);
            delta = atan(params.RPipe/2/params.RPosPipe2);
            xm = ([cos(delta), -sin(delta); sin(delta), cos(delta)]*(p1+p2)'/2)';
            PBC1_1 = p1; % store points for later boundary conditions
            PBC1_2 = p2;
        elseif iPipe == nBars
            mi_drawarc(p2, p1, 180, 20);
            delta = -atan(params.RPipe/2/params.RPosPipe1);
            xm = ([cos(delta), -sin(delta); sin(delta), cos(delta)]*(p1+p2)'/2)';
            PBC2_1 = p1; % store points for later boundary conditions
            PBC2_2 = p2;
        else
            mi_drawarc(p1, p2, 180, 20);
            mi_drawarc(p2, p1, 180, 20);
            xm = (p1+p2)/2;
        end
        mi_addblocklabel(xm);
        mi_selectlabel(xm);
        mi_setblockprop('Air',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
        mi_clearselected();
    end

    mi_getmaterial('M-27 Steel');
    xm = (params.RRI+params.RPosPipe1-params.RPipe)/2*[cos(rotorAngle/2), sin(rotorAngle/2)];
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('M-27 Steel',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    mi_getmaterial('Air');
    xm = (params.RRO+params.RAG)/2*[cos(rotorAngle/2), sin(rotorAngle/2)];
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('Air',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    % BOUNDARY CONDITIONS (for correct meshing)
    mi_addboundprop('Dirichlet',0,0,0,0);
    mi_selectarcsegment(params.RRI*[cos(rotorAngle/2), sin(rotorAngle/2)]);
    mi_setarcsegmentprop(10,'Dirichlet',0,0);
    mi_clearselected

    % Periodic
    mi_addboundprop('AP1',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((point{3}+PBC1_2)/2);
    mi_selectsegment((point{6}+PBC2_2)/2);
    mi_setsegmentprop('AP1',0,1,0,0);
    mi_clearselected
    mi_addboundprop('AP2',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((PBC1_1+PBC1_2)/2);
    mi_selectsegment((PBC2_1+PBC2_2)/2);
    mi_setsegmentprop('AP2',0,1,0,0);
    mi_clearselected
    mi_addboundprop('AP3',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((PBC1_1+point{2})/2);
    mi_selectsegment((PBC2_2+point{5})/2);
    mi_setsegmentprop('AP3',0,1,0,0);
    mi_clearselected

    mi_addboundprop('AP4',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((point{1}+point{2})/2);
    mi_selectsegment((point{4}+point{5})/2);
    mi_setsegmentprop('AP4',0,1,0,0);
    mi_clearselected

    mi_zoomnatural
    
    mi_saveas('TIND_rotor.FEM');

    mi_createmesh();
    mi_analyze(0);

    keyboard
    closefemm();
end