%% Rotor
femm2geo('TIND_rotor.ans', 'TIND_rotor.geo', 'MeshSize', 1e-2, 'RefineArea', [0.185, 0.3], 'RefineRatio', 0.2);
% Rotor bars
nSrf = 0;
bars = [6, 7, 8, 9, 10, 12];
for iBar = 1:10
    nSrf = nSrf + 1;
    barsI = bars+8*(iBar-1);
    barsI(6) = barsI(6)*-1;
    writeSurface('TIND_rotor.geo', 200+nSrf, nSrf, barsI);
end
% Air gap
airGap = [];
for iBar = 1:10
    airGap = [airGap, 91+ (iBar-1)*2];
    airGap = [airGap, [5,12,11] + (iBar-1)*8];
    airGap = [airGap, 92+ (iBar-1)*2];
end
airGap = [airGap, -3, -89, 1];
nSrf = nSrf + 1;
writeSurface('TIND_rotor.geo', 200+nSrf, nSrf, airGap);
% cooling channels
nSrf = nSrf + 1;
writeSurface('TIND_rotor.geo', 200+nSrf, nSrf, [-85,111,112]);
nSrf = nSrf + 1;
writeSurface('TIND_rotor.geo', 200+nSrf, nSrf, [88,149,150]);
for iChannel = 1:9
    nSrf = nSrf + 1;
    writeSurface('TIND_rotor.geo', 200+nSrf, nSrf, [113,114,115,116]+(iChannel-1)*4);
end
% iron body
iron = [];
for iBar = 1:10
    iron = [iron, 91+(iBar-1)*2];
    iron = [iron, [5:11]+(iBar-1)*8];
    iron = [iron, 92+(iBar-1)*2];
end
iron = [iron, 4,88,87,-90,-86,-85,-2];

% iron = [1,2,85,86,90,-87,-88,-4,-3,-89];
nSrf = nSrf + 1;
writeSurface('TIND_rotor.geo', 200+nSrf, nSrf, iron, 212:222);

writePhysicalGroup('TIND_rotor.geo', 'Line', 'Dirichlet', [90]);
writePhysicalGroup('TIND_rotor.geo', 'Line', 'Coupling', [89]);
writePhysicalGroup('TIND_rotor.geo', 'Line', 'AP1', [86,85,2,1]);
writePhysicalGroup('TIND_rotor.geo', 'Line', 'AP2', [87,88,4,3]);

for iBar = 1:10
    writePhysicalGroup('TIND_rotor.geo', 'Surface', ['Bar' +num2str(iBar)], iBar);
end
for iAir = 1:12
    writePhysicalGroup('TIND_rotor.geo', 'Surface', ['Air' +num2str(iAir)], iAir+10);
end
writePhysicalGroup('TIND_rotor.geo', 'Surface', 'Iron', 23);

%% Stator
femm2geo('TIND_stator.ans', 'TIND_stator.geo', 'MeshSize', 1e-2, 'RefineArea', [0.15, 0.195], 'RefineRatio', 0.2);

% Stator slots
nSrf = 0;
slots = [7,15,11,14];
for iSlot = 1:12
    nSrf = nSrf + 1;
    slotI = slots+11*(iSlot-1);
    slotI(4) = slotI(4)*-1;
    writeSurface('TIND_stator.geo', 200+nSrf, nSrf, slotI);
end
slots = [8,9,10,15];
for iSlot = 1:12
    nSrf = nSrf + 1;
    slotI = slots+11*(iSlot-1);
    slotI(4) = slotI(4)*-1;
    writeSurface('TIND_stator.geo', 200+nSrf, nSrf, slotI);
end
% Air gap
airGap = [];
for iSlot = 1:12
    airGap = [airGap, 139+ (iSlot-1)*2];
    airGap = [airGap, [5,6,14,12,13] + (iSlot-1)*11];
    airGap = [airGap, 140+ (iSlot-1)*2];
end
airGap = [airGap, -3, -137, 1];
nSrf = nSrf + 1;
writeSurface('TIND_stator.geo', 200+nSrf, nSrf, airGap);

% cooling channels
for iChannel = 1:15
    nSrf = nSrf + 1;
    writeSurface('TIND_stator.geo', 200+nSrf, nSrf, [163,164,165,166]+(iChannel-1)*4);
end

% iron body
iron = [];
for iSlot = 1:12
    iron = [iron, 139+(iSlot-1)*2];
    iron = [iron, [5:13]+(iSlot-1)*11];
    iron = [iron, 140+(iSlot-1)*2];
end
iron = [iron, 4,-138,-2];

nSrf = nSrf + 1;
writeSurface('TIND_stator.geo', 200+nSrf, nSrf, iron, 226:240);


writePhysicalGroup('TIND_stator.geo', 'Line', 'Dirichlet', [138]);
writePhysicalGroup('TIND_stator.geo', 'Line', 'Coupling', [137]);
writePhysicalGroup('TIND_stator.geo', 'Line', 'AP1', [1,2]);
writePhysicalGroup('TIND_stator.geo', 'Line', 'AP2', [3,4]);

for iSlot = 1:12
    writePhysicalGroup('TIND_stator.geo', 'Surface', ['SlotDown' +num2str(iSlot)], iSlot);
    writePhysicalGroup('TIND_stator.geo', 'Surface', ['SlotUp' +num2str(iSlot)], iSlot+12);
end
for iAir = 25:40
    writePhysicalGroup('TIND_stator.geo', 'Surface', ['Air' +num2str(iAir)], iAir);
end
writePhysicalGroup('TIND_stator.geo', 'Surface', 'Iron', 41);






function writeSurface(name, idLoop, idSurface, lines, idRemove)
    arguments
        name 
        idLoop 
        idSurface 
        lines 
        idRemove = [];
    end
    fid = fopen(name, 'a');
    if ~isempty(idRemove)
        idLoops = [idLoop, idRemove];
    else
        idLoops = idLoop;
    end

    lineStr = join(string(lines), ', ');
    loopStr = join(string(idLoops), ', ');
    
    fprintf(fid, 'Line Loop(%d) = {%s};\n', idLoop, lineStr);
    fprintf(fid, 'Plane Surface(%d) = {%s};\n', idSurface, loopStr);

    fclose(fid);
end

function writePhysicalGroup(fileName, type, name, ids)
    fid = fopen(fileName, 'a');
    
    idStr = join(string(ids), ', ');
    
    fprintf(fid, 'Physical %s("%s")= {%s};\n', type, name, idStr);

    fclose(fid);
end
