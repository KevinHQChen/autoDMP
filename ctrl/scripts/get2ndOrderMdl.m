function [G, y_0, x_0, ns, np, theta0, P_0, mdlStr] = get2ndOrderMdl(no, ni, dt)
%% constants for control model
mdlStr = '2ndOrder';
ns = 2*no;

wn = 1; % Hz
zeta = 1; % critically damped
K = 1;
Ksign = diag(ones(no,1)) - 0.5*(eye(no) ~= 1);

a = wn^2 * dt;
b = 2*zeta*wn * dt;
c = K*wn^2*dt;

A_i = [1, dt;
      -a, 1 - b];
B_i = [zeros(1, ni);
       ones(1, ni)];
C_i = [1, 0];

A_0 = repmat({A_i}, 1, no);
A_0 = blkdiag(A_0{:});          % A = ns x ns   (x, x)
B_0 = zeros(ns, ni);            % B = ns x ni   (x, u)
B_0(2:2:end, :) = c*normalize(Ksign);
C_0 = repmat({C_i}, 1, no);
C_0 = blkdiag(C_0{:});          % C = no x ns   (y, x)
D_0 = zeros(no, ni);            % D = no x ni   (y, u)
K_0 = zeros(ns, no);            % K = ns x no   (y, u)
y0_0 = zeros(no, 1);
x0_0 = zeros(ns, 1);
x0_0(1:2:end) = y0_0;

% double x and y to provide smooth switching between multiple measurements in each channel
x_0 = [x0_0; x0_0];
y_0 = [y0_0; y0_0];
G = idss(blkdiag(A_0, A_0),...  % A = 2*ns x 2*ns   (x, x)
         [B_0; B_0],...         % B = 2*ns x ni     (x, u)
         blkdiag(C_0, C_0),...  % C = 2*no x 2*ns   (y, x)
         [D_0; D_0],...         % D = 2*no x ni     (y, u)
         blkdiag(K_0, K_0),...  % K = 2*ns x 2*no   (x, y)
         x_0, dt);              % x0 = 2*ns x 1     (x, 1)

%% constants for recursive least squares
np = ni + 2; % number of parameters for each output

% initialize paramEst variables
prms = [-a*ones(1, 2*no);
        -b*ones(1, 2*no);
        G.B(2:2:end,:)'];
theta0 = zeros(2*no*np, 1);
for i=1:2*no % stack each col of prms into a single column vector
    theta0( ((i-1)*np + 1):i*np ) = prms(:,i);
end
P_0 = 1*eye(2*no*np);
