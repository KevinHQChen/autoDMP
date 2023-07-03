clc, close all, clear all
rootPath = "~/autoDMP/ctrl/scripts/";
cd(rootPath);
plantConstants;
%% define constants for control

% overall plant model (we can decompose all channel topologies such that each pump has a channel it's directly connected to, and all intermediate channels are just extensions of these direct channels)
A_0 = diag(ones(no,1)-dt);
B_0 = dt*normalize( diag(ones(no,1)) - 0.5*(eye(no) ~= 1) );
C_0 = eye(no);
D_0 = zeros(no, ni);
K_0 = zeros(no);
y0_0 = zeros(no, 1);
x0_0 = y0_0;
% double each output to provide smooth switching between multiple measurements in each channel
G = idss(blkdiag(A_0, A_0),... % A = 2*no x 2*no   (x, x)
         [B_0; B_0],...        % B = 2*no x ni     (x, u)
         blkdiag(C_0, C_0),... % C = 2*no x 2*no   (y, x)
         [D_0; D_0],...        % D = 2*no x ni     (y, u)
         blkdiag(K_0, K_0),... % K = 2*no x 2*no   (x, y)
         [x0_0; x0_0], dt);      % x0 = 2*no x 1     (x, 1)

% initialize normalized MPC weight vectors
uwt0 = zeros(1, ni);
duwt0 = dt*ones(1, ni);
ywt0 = zeros(1, 2*no);

%% define non-virtual buses for AMPC
mdl0 = struct('A', A_0, 'B', B_0, 'C', C_0, 'D', D_0, 'U', u_o, 'Y', y0_0, 'X', x0_0, 'DX', zeros(no, 1));
busInfo = Simulink.Bus.createObject(mdl0);
mdl0_bus = evalin('base', busInfo.busName);
evalin('base', ['clear ' busInfo.busName]);

mdlFull = struct('A', G.A, 'B', G.B, 'C', G.C, 'D', G.D, 'U', u_o, 'Y', y_o, 'X', y_o, 'DX', zeros(2*no, 1));
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
ampc.Model.Nominal.U = u_o;
ampc.Model.Nominal.Y = y_o;
%% specify constraints for MV and MV Rate
for i = 1:ni
    ampc.MV(i).Min = 0;
    ampc.MV(i).Max = u_max(i);
end
%% specify target for MV
for i = 1:ni
    ampc.MV(i).Target = u_o(i);
end
%% specify constraints for OV
for i = 1:2*no
    ampc.OV(i).Min = -y_max(i);
    ampc.OV(i).Max = y_max(i);
end
%% specify overall adjustment factor applied to weights
beta = 0.13534; % maximize robustness in closed-loop performance
%% specify weights
ampc.Weights.MV = uwt0*beta;
ampc.Weights.MVRate = duwt0/beta;
ampc.Weights.OV = ywt0*beta;
ampc.Weights.ECR = 100000;

%% constants for recursive least squares
np = ni+1; % number of parameters for each output

% parameter projection
tau = 10; % memory horizon (in seconds)
forgettingFactor = 1 - dt/tau;

% initialize paramEst variables
prms = [-dt*ones(1, 2*no); G.B'];
theta0 = zeros(2*np*no, 1);
for i=1:2*no % stack each row of prms into a single column vector
    theta0( ((i-1)*np + 1):i*np ) = prms(:,i);
end
P_0 = 1*eye(2*np*no);

%% use custom state estimator implementation
setEstimator(ampc, 'custom');

%% get output disturbance model
s = tf('s');
% 1/s^2 * sign(G.B); % OD model for cross-coupled ramp-like disturbances that account for geometric constraints
% 1/s^2 * diag(ones(numChs,1)); % default OD model for ramp-like disturbances
% tfOD = tf(zeros(2*no));
tfOD = blkdiag(1/s^2 * sign(B_0), 1/s^2 * sign(B_0));
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
    % 1. pre-gen: get in position for droplet generation
    struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', tsl_*5, 'postT', tsl_*5);
    % 2. gen: perform droplet generation
    struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', tsl_*2, 'postT', 0);
    % 3. post-gen: move droplet out of the way for next droplet
    struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', tsl_*5, 'postT', tsl_*5);
    % post-gen cont.: get in position for next droplet generation
    struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', tsl_*2, 'postT', 0);

    struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
    struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);

    struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
    struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);

    struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
    struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);

    struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
    struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);

    struct('r', [0; -Wch/2;  -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; Dneck; -Ld; 0;     0; 0],      'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
    struct('r', [0; 0;       0;  -Wch/2; 0; Ld/2-1], 'preT', 0, 'moveT', ts_*5, 'postT', ts_*5);
    struct('r', [0; 0;       0;  Dplug; 0; Ld/2-1],  'preT', Tpre, 'moveT', ts_*2, 'postT', 0);
             ];

%% trajectory parameters
maxTrajLength = 60/dt; % 60 seconds
tend = 0;
for i=1:length(eventQueue)
    tend = tend + eventQueue(i).preT + eventQueue(i).moveT + eventQueue(i).postT;
end
tend = tend+5;

t = (0:dt:tend)';
r = t;

%% run simulation and analyze results
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
% set_param("T_junction/Selector", 'IndexOptionArray', { 'Select all' }); % capture data from all channels

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
