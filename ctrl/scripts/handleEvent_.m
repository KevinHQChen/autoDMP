function [evDone, waypt, postT] = handleEvent_(event, waypt, postT, y, ymax, numWaypts, dt)
%#codegen
% initialize flags
evDone = false;

% increment waypoint
if waypt < numWaypts
    waypt = waypt + 1;
else
    postT = postT - dt;
end

if postT < 0
    evDone = true;
end
