% 

function params = TIND_stator(params)
    arguments
        % Standard parameters   
        params.Symmetry = 4;
        params.PIPES = 60;              % Number of cooling pipes
        params.SLOTS = 48;              % Number of stator slots
        params.MeshSize = 0;            % Automesh
        params.MeshSizeAG = 1e-3;        % mesh sides at air gap
    
        params.HS = 35.2e-3;            % Height Slot (or depth)
        params.HSO = 4.0e-3;            % Height Slot Opening
        params.RAG = 186.75e-3;         % Radius Air Gap
        params.RPosPipe = 270e-3;       % Radius to center of pipe
        params.RPipe = 7.5e-3;          % Radius of cooling pipes
        params.SP = 1.2e-3;             % Slot Padding
        params.SRI = 187.5e-3;          % Stator Radius Inner
        params.SRO = 300e-3;            % Stator Radius Outer
        params.WS = 10.7e-3;            % Width Slot
        params.WSO = 10.7e-3;           % Width Slot Opening
    end

    % calculated parameters
    a1 = asin((params.WSO/2)/(params.SRI));
    a2 = pi/params.SLOTS;
    a3 = a2-a1;
    slotAngle = 2*pi/params.SLOTS;
    statorAngle = slotAngle*params.SLOTS/params.Symmetry;

    point{1} = params.RAG*[1, 0];
    point{2} = params.SRI*[1, 0]; 
    point{3} = params.SRO*[1, 0];
    point{4} = params.RAG*[cos(statorAngle), sin(statorAngle)];
    point{5} = params.SRI*[cos(statorAngle), sin(statorAngle)];
    point{6} = params.SRO*[cos(statorAngle), sin(statorAngle)];

    point{7} =  [cos(a1)*params.SRI, -params.WSO/2];
    point{8} =  [point{7}(1)+params.HSO, -params.WS/2];
    point{9} =  [point{8}(1), -params.WS/2 - params.SP];
    point{10} =  point{8} + [params.HS/2, 0];
    point{11} =  point{8} + [params.HS, 0];
    point{12} =  point{10} + [0, params.WS];
    point{13} =  point{11} + [0, params.WS];

    point{14} =  [point{9}(1), -point{9}(2)];
    point{15} =  [point{8}(1), -point{8}(2)];
    point{16} =  [point{7}(1), -point{7}(2)];
    point{17} = params.SRI*[cos(slotAngle), sin(slotAngle)];
    HY = params.SRO - point{11}(1);
    % correct unrotated points
    for iPoint = 7:16
        point{iPoint} = ([cos(slotAngle/2), -sin(slotAngle/2); sin(slotAngle/2), cos(slotAngle/2)]*point{iPoint}')';
    end
    nPoints = numel(point);
    % copy and rotate points for all slots
    for iSlot = 1:(params.SLOTS /params.Symmetry)
        for iPoint = 1:nPoints
            rotAngle = slotAngle*iSlot;
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
    mi_drawarc(point{1}, point{4}, rad2deg(statorAngle), 360*params.MeshSizeAG/params.RAG);
    mi_drawarc(point{3}, point{6}, rad2deg(statorAngle), 1);

    mi_drawline(point{1}, point{2});
    mi_drawline(point{2}, point{3});
    mi_drawline(point{4}, point{5});
    mi_drawline(point{5}, point{6});


    % BOUNDARY CONDITIONS 
    mi_addboundprop('Dirichlet',0,0,0,0);
    mi_selectarcsegment(params.SRO*[cos(statorAngle/2), sin(statorAngle/2)]);
    mi_setarcsegmentprop(10,'Dirichlet',0,0);
    mi_clearselected

    mi_addboundprop('AP1',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((point{1}+point{2})/2);
    mi_selectsegment((point{4}+point{5})/2);
    mi_setsegmentprop('AP1',0,1,0,0);
    mi_clearselected
    mi_addboundprop('AP2',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((point{2}+point{3})/2);
    mi_selectsegment((point{5}+point{6})/2);
    mi_setsegmentprop('AP2',0,1,0,0);
    mi_clearselected

    
    mi_getmaterial('10 AWG');

    % create slots
    for iSlot = 1:(params.SLOTS/params.Symmetry)
        mi_drawarc(point{2+(iSlot-1)*nPoints}, point{7+(iSlot-1)*nPoints}, rad2deg(a3), 360*params.MeshSizeAG/params.SRI);
        mi_drawline(point{7+(iSlot-1)*nPoints}, point{9+(iSlot-1)*nPoints});
        mi_drawline(point{9+(iSlot-1)*nPoints}, point{8+(iSlot-1)*nPoints});
        mi_drawline(point{8+(iSlot-1)*nPoints}, point{10+(iSlot-1)*nPoints});
        mi_drawline(point{10+(iSlot-1)*nPoints}, point{11+(iSlot-1)*nPoints});
        mi_drawline(point{11+(iSlot-1)*nPoints}, point{13+(iSlot-1)*nPoints});
        mi_drawline(point{13+(iSlot-1)*nPoints}, point{12+(iSlot-1)*nPoints});
        mi_drawline(point{12+(iSlot-1)*nPoints}, point{15+(iSlot-1)*nPoints});
        mi_drawline(point{15+(iSlot-1)*nPoints}, point{14+(iSlot-1)*nPoints});
        mi_drawline(point{14+(iSlot-1)*nPoints}, point{16+(iSlot-1)*nPoints});
        mi_drawline(point{8+(iSlot-1)*nPoints}, point{15+(iSlot-1)*nPoints});
        mi_drawline(point{10+(iSlot-1)*nPoints}, point{12+(iSlot-1)*nPoints});
        mi_drawarc(point{16+(iSlot-1)*nPoints}, point{17+(iSlot-1)*nPoints}, rad2deg(a3), 360*params.MeshSizeAG/params.SRI);
        % set copper
        xm = (point{8+(iSlot-1)*nPoints} + point{12+(iSlot-1)*nPoints})/2;
        mi_addblocklabel(xm);
        mi_selectlabel(xm);
        mi_setblockprop('10 AWG',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
        mi_clearselected();
        xm = (point{10+(iSlot-1)*nPoints} + point{13+(iSlot-1)*nPoints})/2;
        mi_addblocklabel(xm);
        mi_selectlabel(xm);
        mi_setblockprop('10 AWG',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
        mi_clearselected();
    end

    % create cooling pipes
    for iPipe = 1:(params.PIPES/params.Symmetry)
        rotAngle = -pi/params.PIPES + iPipe * 2*pi/params.PIPES;
        p1 = ([cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)]*[params.RPosPipe + params.RPipe; 0])';
        p2 = ([cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)]*[params.RPosPipe - params.RPipe; 0])';
        mi_drawarc(p1, p2, 180, 20);
        mi_drawarc(p2, p1, 180, 20);
        xm = (p1+p2)/2;
        mi_addblocklabel(xm);
        mi_selectlabel(xm);
        mi_setblockprop('Air',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
        mi_clearselected();
    end

    mi_getmaterial('M-27 Steel');
    xm = (params.SRO+params.RPosPipe+params.RPipe)/2*[cos(statorAngle/2), sin(statorAngle/2)];
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('M-27 Steel',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    mi_getmaterial('Air');
    xm = (params.SRI+params.RAG)/2*[cos(statorAngle/2), sin(statorAngle/2)];
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('Air',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    mi_zoomnatural
    
    mi_saveas('TIND_stator.FEM');
    mi_createmesh();
    mi_analyze(0);

    % keyboard
    closefemm();
end