function [traj, numWaypts] = trajGen_(event, y, ymax, dt, no)
%#codegen
traj = zeros(2*no, 2400);
yDest = zeros(2*no, 1);
directChs = find(event.r); % direct channels
inferredChs = setdiff(1:2*no, directChs); % inferred channels

numWaypts = cast(event.moveT/dt, "uint16");
assert(numWaypts < 1200);

for i=1:length(directChs)
    yDest(directChs(i)) = ymax(directChs(i)).*event.r(directChs(i));
end
for i=1:length(inferredChs) % apply KCL/conservation of charge at junction
    yDest(inferredChs(i)) = -sum(yDest(directChs)) / length(inferredChs);
end

for i=1:size(y,1)
    traj(i, 1:numWaypts) = trapveltraj(...
        [y(i), yDest(i)],...
        numWaypts,...
        PeakVelocity=abs( y(i) - yDest(i) ) / event.moveT);
end
