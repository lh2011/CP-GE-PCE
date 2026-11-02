function [xRot, yRot] = rotdata(xData, yData, rot)
    xRot = xData*cos(rot) - yData*sin(rot);
    yRot = xData*sin(rot) + yData*cos(rot);
end