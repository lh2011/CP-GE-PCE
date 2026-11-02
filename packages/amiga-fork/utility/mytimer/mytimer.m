function [accTime, lastTime] = mytimer(name, status)
    arguments
        name string = ""
        status {mustBeMember(status, ["start", "stop", "print"])} = "print";
    end
    persistent timerStruct;

    if name == "reset"
        timerStruct = struct("Name", {}, "StartTime", {}, "StopTime", {}, "AccumulateTime", {}, "LastTime", {}, "Calls", {});
        return
    end
    
    if isempty(timerStruct)
        pos = 0;
    else
        pos = find(vertcat(timerStruct.Name) == name);
        if isempty(pos) && name ~= ""
            pos = 0;
        end
    end
    
    if pos == 0
        pos = numel(timerStruct) + 1;
        timerStruct(pos).Name = name;
        timerStruct(pos).StartTime = 0;
        timerStruct(pos).StopTime = 0;
        timerStruct(pos).AccumulateTime = 0;
        timerStruct(pos).LastTime = 0;
        timerStruct(pos).Calls = 0;
    end
    
    switch status
        case "start"
            timerStruct(pos).StartTime = tic();
        case "stop"
            timerStruct(pos).StopTime = tic();
            timerStruct(pos).LastTime = toc(timerStruct(pos).StartTime);
            timerStruct(pos).AccumulateTime = timerStruct(pos).AccumulateTime + timerStruct(pos).LastTime;
            timerStruct(pos).Calls = timerStruct(pos).Calls + 1;
        case "print"
            TimerName = vertcat(timerStruct.Name);
            LastTime = vertcat(timerStruct.LastTime);
            AccumulateTime = vertcat(timerStruct.AccumulateTime);
            AverageTime = vertcat(timerStruct.AccumulateTime)./vertcat(timerStruct.Calls);
            Calls = vertcat(timerStruct.Calls);
            timerTable = table(TimerName, LastTime, AccumulateTime, AverageTime, Calls);
            disp(timerTable);
    end

    if nargout
        lastTime = timerStruct(pos).LastTime;
        accTime = timerStruct(pos).AccumulateTime;
    end
end