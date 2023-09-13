function [A, B, C, D, Q, R, N] = stateEst_(Ap, Bp, Cp, Dp, Aod, Bod, Cod, Dod, Dn, no, ni, ns_)
%#codegen
nsp = ns_;             % n_plant_states
nsod = size(Aod,1);     % n_od_states
ns = nsp + nsod;        % n_states = n_plant_states + n_od_states

A = zeros(ns);      % n_states x n_states
B = zeros(ns,ni);   % n_states  x n_inputs
C = zeros(no,ns);   % n_outputs x n_states
D = zeros(no,ni);   % n_outputs x n_inputs
Q = zeros(ns,ns);   % n_states  x n_states
G = eye(ns);        % n_states  x n_states
R = zeros(no,no);   % n_outputs x n_outputs
N = zeros(ns,no);   % n_states  x n_outputs
H = zeros(no,ns);   % n_outputs x n_states

% combine plant and output disturbance model
% (force the outputs to fit in preallocated memory)
A(1:ns, 1:ns) = blkdiag(Ap, Aod);
B(1:nsp, 1:ni) = Bp;
C(1:no, 1:ns) = [Cp Cod];
D(1:no, 1:ni) = Dp;

B_est = zeros(ns, ni + no + no);
B_est(1:ns, 1:ni+no) = blkdiag(Bp, Bod);
D_est = [Dp Dod Dn];
Q = B_est * B_est';
R = D_est * D_est';
N = B_est * D_est';

% [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
% [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);

% xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
% yhat = C*xhat + D*u;
