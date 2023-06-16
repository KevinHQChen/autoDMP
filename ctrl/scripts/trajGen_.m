function [traj, numWaypts] = trajGen_(event, y, ymax, dt)
%#codegen
traj = zeros(size(y,1), 2400);
yDest = zeros(size(y));
chs = find(event.r); % directly measured channels
invChs = setdiff(1:size(y,1), chs); % indirectly measured channels

numWaypts = cast(event.moveT/dt, "uint16");
assert(numWaypts < 1200);

% generate trajectory for directly measured channels
for i=1:length(chs)
    yDest(chs(i)) = ymax(chs(i)).*event.r(chs(i));
    traj(chs(i), 1:numWaypts) = trapveltraj(...
        [y(chs(i)), yDest(chs(i))],...
        numWaypts,...
        PeakVelocity=abs( y(chs(i)) - yDest(chs(i)) ) / event.moveT);
end

% generate trajectory for indirectly measured channels
for i=1:length(invChs)
    yDest(invChs(i)) = -sum(yDest(chs)) / length(invChs); % apply KCL to estimate unmeasurable channels
    traj(invChs(i), 1:numWaypts) = trapveltraj(...
        [y(invChs(i)), yDest(invChs(i))],...
        numWaypts,...
        PeakVelocity=abs( y(invChs(i)) - yDest(invChs(i)) ) / event.moveT);
end
