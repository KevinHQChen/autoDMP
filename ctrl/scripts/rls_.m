function [dtheta, dP, L] = rls_(theta, phi, epsil, EN, p_, dPmod_, lambda, P)
%#codegen
no = size(phi,2);
ni = no;
np = ni + 1;
dtheta = zeros(np*no, 1);
dP = zeros(np*no,np*no);
L = zeros(np*no, no);

L = P*phi*inv(lambda*eye(no) + phi'*P*phi);

dtheta(1:np*no, 1) = L*epsil;
dP(1:np*no, 1:np*no) = L*phi'*P;

% parameter projection
for i = 1:np*no
    if mod(i-1, np) == 0 % a_i requires lower and upper bounds
        if ~( (theta(i) + dtheta(i) > p_) | ... % lower bound
              ( (theta(i) + dtheta(i) == p_) & (dtheta(i) >= 0) ) ) ...
         & ~( (theta(i) + dtheta(i) < p_+0.5) | ... % upper bound
              ( (theta(i) + dtheta(i) == p_+0.5) & (dtheta(i) <= 0) ) )
            dtheta(i) = 0;
            dP(i,i) = 0;
            dP(i,i) = 0;
            % dP(i, i:i+np-1) = 0;
            % dP(i:i+np-1, i) = 0;
        end
    else % b_ij only requires lower bound
        if ~( (theta(i) + dtheta(i) > p_) || ...
              ( (theta(i) + dtheta(i) == p_) && (dtheta(i) >= 0) ) )
            dtheta(i) = 0;
            dP(i,i) = 0;
            dP(i,i) = 0;
        end
    end
end

for i = 1:no
    if ~EN(i)
        dtheta( ((i-1)*np + 1):i*np ) = 0;
        dP( ((i-1)*np + 1):i*np,: ) = 0;
        dP( :,((i-1)*np + 1):i*np ) = 0;
    end
end

% if ~( all(theta+dtheta > p_) || ...
%       (all(theta+dtheta >= p_) && all(b_(find(a_)))) )
%     dtheta = zeros(numChs,numChs);
%     dP = zeros(numChs,numChs);
%     % dP = - dPmod_*eye(numChs); % covariance "kick"
% end

% covariance resetting (doesn't work well)
% if norm(P - dP) < 5e-2
%     dP = - ( P + 5e-1*eye(3) );
% end
% if min(svd(P - dP)) < 1e-3
%     dP = - ( P + 1*eye(3) );
% end
