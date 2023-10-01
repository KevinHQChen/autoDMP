clc, close all, clear all
rootPath = "~/autoDMP/ctrl/scripts/";
cd(rootPath);
plantConstants;

%% define overall plant model
[G, y_0, x_0, ns, np, theta0, P_0, mdlStr] = get1stOrderMdl(no, ni, dt);
% [G, y_0, x_0, ns, np, theta0, P_0, mdlStr] = get2ndOrderMdl(no, ni, dt);

keySet = {'integrator', '1stOrder', '2ndOrder'};
valSet = {0, 1, 2};
mdlMap = containers.Map(keySet, valSet);
mdlNum = mdlMap(mdlStr);

%% define constants for control
G0 = load('G0').G;
G1 = load('G1').G;
G2 = load('G2').G;
Gg = load('GgFull').G;
% G0 = load('simSysID/G0').G;
% G1 = load('simSysID/G1').G;
% G2 = load('simSysID/G2').G;

%% define non-virtual buses for AMPC
mdlFull = struct('A', G.A, 'B', G.B, 'C', G.C, 'D', G.D, 'U', u_0, 'Y', y_0, 'X', x_0, 'DX', zeros(2*ns, 1));
busInfo = Simulink.Bus.createObject(mdlFull);
mdl_bus = evalin('base', busInfo.busName);
evalin('base', ['clear ' busInfo.busName]);

%% create MPC controller object with sample time
ampc = mpc(G, dt);
%% specify prediction horizon
ampc.PredictionHorizon = 20;
%% specify control horizon
ampc.ControlHorizon = 1;
%% specify nominal values for inputs and outputs
ampc.Model.Nominal.U = u_0;
ampc.Model.Nominal.Y = y_0;
%% specify constraints for MV and MV Rate
for i = 1:ni
    ampc.MV(i).Min = 0;
    ampc.MV(i).Max = u_max(i);
end
%% specify target for MV
for i = 1:ni
    ampc.MV(i).Target = u_0(i);
end
%% specify constraints for OV
for i = 1:2*no
    ampc.OV(i).Min = -y_max(i);
    ampc.OV(i).Max = y_max(i);
end
%% specify overall adjustment factor applied to weights
beta = 0.13534; % maximize robustness in closed-loop performance
%% initialize normalized MPC weight vectors
uwt0 = zeros(1, ni);
duwt0 = dt*ones(1, ni);
ywt0 = zeros(1, 2*no);
%% specify weights
ampc.Weights.MV = uwt0*beta;
ampc.Weights.MVRate = duwt0/beta;
ampc.Weights.OV = ywt0*beta;
ampc.Weights.ECR = 100000;

%% create MPC controller object with sample time
% add integral action
A = [G0.A, zeros(size(G0.A, 1), 1);
     dt, 1];
B = [G0.B;
     zeros(1, 3)];
C = [G0.C, zeros(size(G0.C, 1), 1);
     zeros(size(G0.C, 2), size(G0.A, 1)), 1];
D = [G0.D;
     zeros(1, 3)];
G0i = ss(A, B, C, D, dt);
mpc1 = mpc(G0i, dt);
%% specify prediction horizon
mpc1.PredictionHorizon = 20;
%% specify control horizon
mpc1.ControlHorizon = 1;
%% specify nominal values for inputs and outputs
mpc1.Model.Nominal.U = u_0;
mpc1.Model.Nominal.Y = [0;0;0;0];
%% specify constraints for MV and MV Rate
for i = 1:ni
    mpc1.MV(i).Min = 0;
    mpc1.MV(i).Max = u_max(i);
end
%% specify constraints for OV
for i = 1:no
    mpc1.OV(i).Min = -y_max(i);
    mpc1.OV(i).Max = y_max(i);
end
mpc1.OV(4).Min = -y_max(1)*dt;
mpc1.OV(4).Max = y_max(1)*dt;
%% specify overall adjustment factor applied to weights
beta = 0.13534; % maximize robustness in closed-loop performance
%% specify weights
mpc1.Weights.MV = uwt0*beta;
mpc1.Weights.MVRate = duwt0/beta;
mpc1.Weights.OV = [1, 0, 0, 1]*beta;
mpc1.Weights.ECR = 100000;

%% use custom state estimator implementation
s = tf('s');
tfOD = [1/s^2 * (eye(3) - (eye(no) ~= 1)); zeros(1, 3)];
setoutdist(mpc1, 'model', tfOD);
God1 = getoutdist(mpc1);
Aod1 = God1.A;
Bod1 = God1.B;
Cod1 = God1.C;
Dod1 = God1.D;
setEstimator(mpc1, 'custom');

%% create MPC controller object with sample time
% add integral action
A = [G1.A, zeros(size(G1.A, 1), 2);
     dt*eye(2), eye(2)];
B = [G1.B;
     zeros(2, 3)];
C = [G1.C, zeros(size(G1.C, 1), 2);
     zeros(size(G1.C, 2), size(G1.A, 1)), eye(2)];
D = [G1.D;
     zeros(2, 3)];
G1i = ss(A, B, C, D, dt);
mpc2 = mpc(G1i, dt);
%% specify prediction horizon
mpc2.PredictionHorizon = 20;
%% specify control horizon
mpc2.ControlHorizon = 1;
%% specify nominal values for inputs and outputs
mpc2.Model.Nominal.U = u_0;
mpc2.Model.Nominal.Y = [0;0;0;0;0];
%% specify constraints for MV and MV Rate
for i = 1:ni
    mpc2.MV(i).Min = 0;
    mpc2.MV(i).Max = u_max(i);
end
%% specify constraints for OV
for i = 1:no
    mpc2.OV(i).Min = -y_max(i);
    mpc2.OV(i).Max = y_max(i);
end
mpc2.OV(4).Min = -y_max(2)*dt;
mpc2.OV(4).Max = y_max(2)*dt;
mpc2.OV(5).Min = -y_max(3)*dt;
mpc2.OV(5).Max = y_max(3)*dt;
%% specify overall adjustment factor applied to weights
beta = 0.13534; % maximize robustness in closed-loop performance
%% specify weights
mpc2.Weights.MV = uwt0*beta;
mpc2.Weights.MVRate = duwt0/beta;
mpc2.Weights.OV = [0, 1, 1, 0, 0]*beta;
mpc2.Weights.ECR = 100000;

%% use custom state estimator implementation
tfOD = [1/s^2 * (eye(3) - (eye(no) ~= 1)); zeros(2, 3)];
setoutdist(mpc2, 'model', tfOD);
God2 = getoutdist(mpc2);
Aod2 = God2.A;
Bod2 = God2.B;
Cod2 = God2.C;
Dod2 = God2.D;
setEstimator(mpc2, 'custom');

%% create MPC controller object with sample time
% add integral action
A = [G2.A, zeros(size(G2.A, 1), 2);
     dt*eye(2), eye(2)];
B = [G2.B;
     zeros(2, 3)];
C = [G2.C, zeros(size(G2.C, 1), 2);
     zeros(size(G2.C, 2), size(G2.A, 1)), eye(2)];
D = [G2.D;
     zeros(2, 3)];
G2i = ss(A, B, C, D, dt);
mpc3 = mpc(G2i, dt);
%% specify prediction horizon
mpc3.PredictionHorizon = 20;
%% specify control horizon
mpc3.ControlHorizon = 1;
%% specify nominal values for inputs and outputs
mpc3.Model.Nominal.U = u_0;
mpc3.Model.Nominal.Y = [0;0;0;0;0];
%% specify constraints for MV and MV Rate
for i = 1:ni
    mpc3.MV(i).Min = 0;
    mpc3.MV(i).Max = u_max(i);
end
%% specify constraints for OV
for i = 1:no
    mpc3.OV(i).Min = -y_max(i);
    mpc3.OV(i).Max = y_max(i);
end
% mpc3.OV(4).Min = -y_max(1)*dt;
% mpc3.OV(4).Max = y_max(1)*dt;
% mpc3.OV(5).Min = -y_max(2)*dt;
% mpc3.OV(5).Max = y_max(2)*dt;
%% specify overall adjustment factor applied to weights
beta = 0.13534; % maximize robustness in closed-loop performance
%% specify weights
mpc3.Weights.MV = uwt0*beta;
mpc3.Weights.MVRate = duwt0/beta;
mpc3.Weights.OV = [1, 1, 0, 1, 1]*beta;
mpc3.Weights.ECR = 100000;

%% use custom state estimator implementation
tfOD = [1/s^2 * (eye(3) - (eye(no) ~= 1)); zeros(2, 3)];
setoutdist(mpc3, 'model', tfOD);
God3 = getoutdist(mpc3);
Aod3 = God3.A;
Bod3 = God3.B;
Cod3 = God3.C;
Dod3 = God3.D;
setEstimator(mpc3, 'custom');

%% create MPC controller object with sample time
% add integral action
A = [G2.A, zeros(size(G2.A, 1), 2);
     dt*eye(2, size(G2.A, 1)), eye(2)];
B = [G2.B;
     zeros(2, 3)];
C = [G2.C, zeros(size(G2.C, 1), 2);
     zeros(2, size(G2.A, 1)), eye(2)];
D = [G2.D;
     zeros(2, 3)];
G2i = ss(A, B, C, D, dt);
mpcg = mpc(G2i, dt);
%% specify prediction horizon
mpcg.PredictionHorizon = 20;
%% specify control horizon
mpcg.ControlHorizon = 1;
%% specify nominal values for inputs and outputs
mpcg.Model.Nominal.U = u_0;
mpcg.Model.Nominal.Y = [0;0;0;0;0];
%% specify constraints for MV and MV Rate
for i = 1:ni
    mpcg.MV(i).Min = 0;
    mpcg.MV(i).Max = u_max(i);
end
%% specify constraints for OV
for i = 1:no
    mpcg.OV(i).Min = -y_max(i);
    mpcg.OV(i).Max = y_max(i);
end
mpcg.OV(4).Min = -y_max(1)*dt;
mpcg.OV(4).Max = y_max(1)*dt;
mpcg.OV(5).Min = -y_max(2)*dt;
mpcg.OV(5).Max = y_max(2)*dt;
%% specify overall adjustment factor applied to weights
beta = 0.13534; % maximize robustness in closed-loop performance
%% specify weights
mpcg.Weights.MV = uwt0*beta;
mpcg.Weights.MVRate = duwt0/beta;
mpcg.Weights.OV = [1, 1, 0, 0, 0]*beta;
mpcg.Weights.ECR = 100000;

%% use custom state estimator implementation
% tfOD = [1/s^2 * (eye(3) - (eye(no) ~= 1)); zeros(2, 3)];
% setoutdist(mpcg, 'model', tfOD);
% God3 = getoutdist(mpcg);
% Aod3 = God3.A;
% Bod3 = God3.B;
% Cod3 = God3.C;
% Dod3 = God3.D;
% setEstimator(mpcg, 'custom');

mdl0 = struct('A', G0.A, 'B', G0.B, 'C', G0.C, 'D', G0.D, 'U', u_0, 'Y', [0;0;0], 'X', 0, 'DX', 0);
busInfo = Simulink.Bus.createObject(mdl0);
mdl0_bus = evalin('base', busInfo.busName);
evalin('base', ['clear ' busInfo.busName]);

mdl1 = struct('A', G1.A, 'B', G1.B, 'C', G1.C, 'D', G1.D, 'U', u_0, 'Y', [0;0;0], 'X', [0;0], 'DX', [0;0]);
busInfo = Simulink.Bus.createObject(mdl1);
mdl1_bus = evalin('base', busInfo.busName);
evalin('base', ['clear ' busInfo.busName]);

mdl2 = struct('A', G2.A, 'B', G2.B, 'C', G2.C, 'D', G2.D, 'U', u_0, 'Y', [0;0;0], 'X', [0;0], 'DX', [0;0]);
busInfo = Simulink.Bus.createObject(mdl2);
mdl2_bus = evalin('base', busInfo.busName);
evalin('base', ['clear ' busInfo.busName]);

%% get noise covariances
Gmn1 = ss(eye(3), Ts=dt); % default measurement noise model
Dmn1 = Gmn1.D;

%% constants for recursive least squares
% parameter projection
tau = 10; % memory horizon (in seconds)
forgettingFactor = 1 - dt/tau;

%% use custom state estimator implementation
setEstimator(ampc, 'custom');

%% get output disturbance model
s = tf('s');
% 1/s^2 * sign(G.B); % OD model for cross-coupled ramp-like disturbances that account for geometric constraints
% 1/s^2 * diag(ones(numChs,1)); % default OD model for ramp-like disturbances
% tfOD = tf(zeros(2*no));
B_0 = G.B(1:ns, :);
if mdlNum == 0 || mdlNum == 1
    tfOD = blkdiag(1/s^2 * sign(B_0), 1/s^2 * sign(B_0));
elseif mdlNum == 2
    tfOD = blkdiag( 1/s^2 * sign(B_0(2:2:end,:)), 1/s^2 * sign(B_0(2:2:end,:)) );
end
setoutdist(ampc, 'model', tfOD);
God = getoutdist(ampc);
Aod = God.A;
Bod = God.B;
Cod = God.C;
Dod = God.D;

nsod = size(God.A,1);
Aod0 = God.A(1:nsod/2, 1:nsod/2);
Bod0 = God.B(1:nsod/2, 1:ni);
Cod0 = God.C(1:no, 1:nsod/2);
Dod0 = God.D(1:no, 1:ni);

%% get noise covariances
Dmn0 = eye(no);
Gmn = ss(blkdiag(Dmn0, Dmn0), Ts=dt); % default measurement noise model
Dmn = Gmn.D;
B_est = [blkdiag(G.B, God.B), zeros(size(G.B,1) + size(God.B,1), size(G.C,1))];
D_est = [G.D God.D Gmn.D];
Q = B_est * B_est';
R = D_est * D_est';
N = B_est * D_est';

%% lowpass filter for paramEst error
lpf = designfilt('lowpassfir','SampleRate',samplingRate,'PassbandFrequency',5, ...
'StopbandFrequency',7.5,'PassbandRipple',0.1, ...
'StopbandAttenuation',60,'DesignMethod','kaiserwin');
[lpfNum, lpfDen] = tf(lpf);

%% define event bus
nullEv = struct('r', [0; 0; 0; 0; 0; 0], 'preT', 0, 'moveT', 0, 'postT', 0);
busInfo = Simulink.Bus.createObject(nullEv);
event_bus = evalin('base', busInfo.busName);
evalin('base', ['clear ' busInfo.busName]);

%% define event sequence
Ld = 0.25; % droplet length (% of channel length)
Wch = 0.02; % channel width (% of channel length)
Dneck = Wch; % required displacement into each unmeasured channel for droplet split/gen (% of channel length)
Tpre = 0.01; % default hold time before move
tsl_ = 2; % learning timescale
ts_ = 2; % validation timescale
Dplug = Wch/2;

eventQueue = [...
    struct('r', [-0.5;  -0.5;  0; 0; 0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
    struct('r', [-0.25; -0.75; 0; 0; 0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
    struct('r', [-0.5;  -0.75; 0; 0; 0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
    struct('r', [-0.5; -0.5;   0; 0; 0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
             ]

% eventQueue = [...
%     struct('r', [0; 0;  0; -0.5;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     struct('r', [0; 0;  0; -Wch/2;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     struct('r', [0; 0;  0; Wch/2;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     struct('r', [0; 0;  0; -0.5;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     struct('r', [0; 0;  0; -0.75;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     struct('r', [0; 0;  0; -0.5;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     % 3. post-gen: move droplet out of the way for next droplet
%     struct('r', [0; 0;       0;  -Wch/2; 0; 0], 'preT', 0, 'moveT', tsl_*5, 'postT', tsl_*5);
%     % post-gen cont.: get in position for next droplet generation
%     struct('r', [0; 0;       0;  Dplug; 0; 0],  'preT', Tpre, 'moveT', tsl_*2, 'postT', 0);
%     % 1. pre-gen: get in position for droplet generation
%     struct('r', [0; 0;  0; 0;     -Wch/2; -Ld],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     % % 1. pre-gen: get in position for droplet generation
%     % struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%              ]

% eventQueue = [...
%     % 1. pre-gen: get in position for droplet generation
%     struct('r', [0; -200/612;  -84/612; 0;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     % 2. gen: perform droplet generation
%     struct('r', [0; -200/612; 84/612; 0;     0; 0],      'preT', Tpre, 'moveT', tsl_*2, 'postT', 0);
%     % 3. post-gen: move droplet out of the way for next droplet
%     struct('r', [0; 0;       0;  -84/612; -406/612; 0], 'preT', 0, 'moveT', tsl_*5, 'postT', tsl_*5);
%              ];

% eventQueue = [...
%     % 1. pre-gen: get in position for droplet generation
%     struct('r', [0; -Ld;  -Wch/2; 0;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     % 2. gen: perform droplet generation
%     struct('r', [0; -Ld; Dneck; 0;     0; 0],      'preT', Tpre, 'moveT', tsl_*2, 'postT', 0);
%     % 3. post-gen: move droplet out of the way for next droplet
%     struct('r', [0; 0;       0;  -Wch/2; Ld/2-1; 0], 'preT', 0, 'moveT', tsl_*5, 'postT', tsl_*5);
%     % post-gen cont.: get in position for next droplet generation
%     struct('r', [0; 0;       0;  Dplug; Ld/2-1; 0],  'preT', Tpre, 'moveT', tsl_*2, 'postT', 0);

%     struct('r', [0; -Ld;  -Wch/2; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; -Ld; Dneck; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
%     struct('r', [0; 0;       0;  -Wch/2; Ld/2-1; 0], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; 0;       0;  Dplug; Ld/2-1; 0],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
%              ];

% eventQueue = [...
%     % 1. pre-gen: get in position for droplet generation
%     struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
%     % 2. gen: perform droplet generation
%     struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', tsl_*2, 'postT', 0);
%     % 3. post-gen: move droplet out of the way for next droplet
%     struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', tsl_*5, 'postT', tsl_*5);
%     % post-gen cont.: get in position for next droplet generation
%     struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', tsl_*2, 'postT', 0);

%     struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
%     struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);

%     struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
%     struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);

%     struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
%     struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);

%     struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
%     struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);

%     struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
%     struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
%     struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
%              ];

%% trajectory parameters
maxTrajLength = 60/dt; % 60 seconds
tend = 0;
for i=1:length(eventQueue)
    tend = tend + eventQueue(i).preT + eventQueue(i).moveT + eventQueue(i).postT;
end
tend = tend+5;

t = (0:dt:tend)';
r = t;

%% run simulation
sim = 'T_junction_mpc';
open(sim);
commandwindow;

% input data setup
% WARNING for multidimensional inputs we have to manually set inport dimensions
% set_param(sim, 'LoadExternalInput', 'on',...
%                'ExternalInput', '[t, r]');

% output data setup
% WARNING for discrete dynamical systems we need to set sampling rates for the 1/z
% blocks manually (to 40Hz in this case) or we'll end up updating it at the
% solver stepSize even if OutputTimes was set properly
set_param(sim, 'ReturnWorkspaceOutputs', 'on',...
               'ReturnWorkspaceOutputsName', 'out',...
               'SignalLogging', 'on',...
               'SignalLoggingName', 'logsout',...
               'SaveFormat', 'Dataset',...
               'OutputOption', 'SpecifiedOutputTimes',...
               'OutputTimes', 't');

% solver setup
% fixed solver step size mitigates consecutive zero crossing error (likely due to abs values)
% and mitigates extremely small time steps due to certain blocks (e.g. random number block with very short sampling periods)
stepSize = dt; % 1e-3;
set_param(sim, 'StopTime', num2str(t(end)),...
               'SolverType', 'Fixed-step',...
               'Solver', 'FixedStepAuto',...
               'FixedStep', num2str(stepSize));
set_param(sim, 'StopTime', num2str(t(end)));

% run simulation and return logged data to workspace
disp('running simulation...');
set_param(sim, 'SimulationCommand', 'start');
set_param(sim, 'SimulationCommand', 'WriteDataLogs');
simStatus = 'run';

while ~strcmp(simStatus, 'stopped')
    simStatus = get_param(sim, 'SimulationStatus');
    pause(0.01);     % pass thread execution to simulink (https://www.mathworks.com/matlabcentral/answers/290503-pause-matlab-to-wait-for-simulink-in-r2016a)
end
disp('stopped simulation.');
signals = getSignals(sim, out);
disp('signals captured:');
keys(signals)     % list all captured signals

%% plot results
t = signals('t');
y = squeeze(signals('y'));
ywt = squeeze(signals('ywt'));
r = squeeze(signals('r'));
theta = squeeze(out.logsout.getElement(3).Values.Data);

figure
plot(t, y([2 3 4 5],:))
grid
legend('y2', 'y3', 'y4', 'y5');
xlabel('Time [s]')
ylabel('Position [px]')
set(gcf, 'Position', [100, 100, 4*300, 300]);
set(gcf, 'Color', 'w');
% export_fig /home/khqc/thesis/assets/dropgenval.png

figure
plot(t, ywt([2 3 4 5],:))
grid
legend('ywt2', 'ywt3', 'ywt4', 'ywt5');
xlabel('Time [s]')
ylabel('Output weight')
set(gcf, 'Position', [100, 100, 4*300, 300]);
set(gcf, 'Color', 'w');
% export_fig /home/khqc/thesis/assets/dropgenywt.png

figure
plot(t, theta)
grid
xlabel('Time [s]')
ylabel('Parameter Value')
set(gcf, 'Position', [100, 100, 4*300, 300]);
set(gcf, 'Color', 'w');
% export_fig /home/khqc/thesis/assets/dropgentheta.png

% figure
% signals('t') = t; % when using fixed-step solvers
% plot(signals('t'), signals('r'))
% hold on
% plot(signals('t'), signals('y'))
% hold on
% plot(signals('t'), signals('u'))
% grid
