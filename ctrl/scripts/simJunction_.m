function RST = simJunction_(y_, no, currTraj, currEv, ywt, beta)
persistent yPrev rPrev currTrajPrev rstCountDown;
y = [y_(1:2:end); y_(2:2:end)];
RST = logical(zeros(2*no, 1));

% initialize variables
if isempty(yPrev)
    yPrev = y;
end
if isempty(rPrev)
    rPrev = currEv.r;
end
if isempty(currTrajPrev)
    currTrajPrev = currTraj;
end
if isempty(rstCountDown) || any(rPrev ~= currEv.r)
    rstCountDown = 2; % reset hysteresis counter on init and on event transition
end

% if primary measurements are controlled and any zero-crossing occurs,
% reset secondary measurements
if any(currEv.r(1:no)) && any(sign(currTraj(1:no)) ~= sign(currTrajPrev(1:no)))
    RST(2:2:end) = true;
else
    RST(2:2:end) = false;
end
% if secondary measurements are controlled and any zero-crossing occurs,
% reset primary measurements
if any(currEv.r(no+1:2*no)) && any(sign(currTraj(no+1:2*no)) ~= sign(currTrajPrev(no+1:2*no)))
    RST(1:2:end) = true;
else
    RST(1:2:end) = false;
end
% if there is contention in control btwn primary and secondary, don't reset any measurement
if any(ywt < 0.995*beta & ywt > 0.005*beta)
    RST = logical(zeros(2*no, 1));
end
yPrev = y_;
rPrev = currEv.r;
currTrajPrev = currTraj;

% hysteresis
if rstCountDown > 0
    RST = logical(zeros(2*no, 1));
    rstCountDown = rstCountDown - 1;
end
if any(RST)
    rstCountDown = 2;
end
