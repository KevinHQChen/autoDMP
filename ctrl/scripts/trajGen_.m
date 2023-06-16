function [traj, numWaypts] = trajGen_(event, y, ymax, dt, no)
%#codegen
traj = zeros(2*no, 2400);
yDest = zeros(2*no, 1);
chs = find(event.r); % directly measured channels
invChs = setdiff(1:2*no, chs); % indirectly measured channels

numWaypts = cast(event.moveT/dt, "uint16");
assert(numWaypts < 1200);

for i=1:length(chs) % directly measured channels
    yDest(chs(i)) = ymax(chs(i)).*event.r(chs(i));
end
for i=1:length(invChs) % indirectly measured channels (by applying KCL/conservation of charge)
    yDest(invChs(i)) = -sum(yDest(chs)) / length(invChs);
end

for i=1:size(y,1)
    traj(i, 1:numWaypts) = trapveltraj(...
        [y(i), yDest(i)],...
        numWaypts,...
        PeakVelocity=abs( y(i) - yDest(i) ) / event.moveT);
end
