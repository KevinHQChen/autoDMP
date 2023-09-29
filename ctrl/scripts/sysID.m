clc, close all, clear all
rootPath = "~/autoDMP/ctrl/scripts/";
%% global variables
samplingRate = 40;
clockPeriod = 4;
dt = 1/samplingRate;

trainPath = "~/thesis/data/train.csv";
% trainPath = "~/autoDMP/ctrl/scripts/simSysID/train/ctrlDataQueue.txt";
trainStart = 10/dt; % state0
% trainEnd = clockPeriod*2^10-3-1-1; % 102.3/dt
trainEnd = clockPeriod*2^9-3-1-1; % - 12.3/dt; % 102.3/dt (for state1)

valPath = "~/thesis/data/val.csv";
% valPath = "~/autoDMP/ctrl/scripts/simSysID/val/ctrlDataQueue.txt";
valStart = 10/dt;
% valEnd = clockPeriod*2^10-3-1-1; % 102.3/dt
valEnd = clockPeriod*2^9-3-1-1; % - 12.3/dt; % 102.3/dt (for state1)

col = dictionary(["t", "y0", "y1", "y2", "u0", "u1", "u2"], 1:7);

idMdlName = 'G2'; % [G0 | G1 | G2]
inputs = col(["u0" "u1" "u2"]);
inputNames = {'Pump1', 'Pump2', 'Pump3'};

if strcmp(idMdlName, 'G0')
    outputs = col(["y0"]); % ["y0" "y1" "y2"]
    outputNames = {'Ch1Water'}; % {'Ch1Water', 'Ch2Oil', 'Ch3Oil'};
elseif strcmp(idMdlName, 'G1')
    outputs = col(["y1" "y2"]);
    outputNames = {'Ch2Oil', 'Ch3Oil'};
elseif strcmp(idMdlName, 'G2')
    outputs = col(["y0" "y1"]); % ["y0" "y1" "y2"]
    outputNames = {'Ch1Water', 'Ch2Oil'}; % {'Ch1Water', 'Ch2Oil', 'Ch3Oil'};
end
no = length(outputs);

%% load data
% train
data_train = readmatrix(trainPath);
data_train = data_train(1:end,:);
t_train = (0:dt:(size(data_train,1) - 1)*dt)'; % assume uniform sampling
u_train = data_train(:, inputs);
y_train = data_train(:, outputs);
% u_train_ts = timeseries(u_train, t_train);
% y_train_ts = timeseries(y_train, t_train);

% val
data_val = readmatrix(valPath);
data_val = data_val(1:end,:);
t_val = (0:dt:(size(data_val,1) - 1)*dt)'; % assume uniform sampling
u_val = data_val(:, inputs);
y_val = data_val(:, outputs);
% u_val_ts = timeseries(u_val, t_val);
% y_val_ts = timeseries(y_val, t_val);

figure
plot(t_train, y_train, t_val, y_val);
legend([outputNames, outputNames]);
print(rootPath + "plots/raw", "-dpng");
rootPath + "plots/raw" + ".png";

%% preprocess data
% remove outliers
% y_train = filloutliers(y_train, 'linear');
% y_val = filloutliers(y_val, 'linear');

sys_train = iddata(y_train(trainStart:trainEnd,:), u_train(trainStart:trainEnd,:), dt);
sys_val = iddata(y_val(valStart:valEnd,:), u_val(valStart:valEnd,:), dt);
sys_train.InputName  = inputNames;
sys_val.InputName  = inputNames;
sys_train.OutputName = outputNames;
sys_val.OutputName = outputNames;

sys_traind = detrend(sys_train, 0);
% sys_traindf = idfilt(sys_traind,5,0.049975); % 0.1 - 6.28 rad/s
sys_vald = detrend(sys_val, 0);
% sys_valdf = idfilt(sys_vald,5,0.049975); % 0.1 - 6.28 rad/s

sys_traindf = fft(sys_traind);
sys_valdf = fft(sys_vald);

figure
plot(sys_traind)
print(rootPath + "plots/trainIDData", "-dpng");
rootPath + "plots/trainIDData" + ".png";

figure
plot(sys_vald)
print(rootPath + "plots/valIDData", "-dpng");
rootPath + "plots/valIDData" + ".png";

%% process estimation - 1 pole
Opt = procestOptions;
Opt.WeightingFilter = [0 31.4159];
proc_ = procest(sys_traindf,{'P1', 'P1', 'P1'} , Opt);

%% process estimation - 2 pole (underdamped)
Opt = procestOptions;
Opt.WeightingFilter = [0 31.4159];
proc_ = procest(sys_traindf,{'P2U', 'P2U', 'P2U'} , Opt);

%% process estimation - 2 pole (underdamped)
Opt = procestOptions;
proc_ = procest(sys_traind,{'P2U', 'P2U', 'P2U'} , Opt);

%% process estimation - 2 pole, (underdamped), 1 integrator, 1 zero
Opt = procestOptions;
Opt.WeightingFilter = [0 31.4159];
proc_ = procest(sys_traindf,{'P2IZU', 'P2IZU', 'P2IZU'} , Opt);

%% convert to state-space
proc_est_ = minreal(ss(proc_));

if rank(ctrb(proc_est_)) < size(proc_est_.A, 1)
    disp('Plant is not controllable');
else
    disp('Plant is controllable');
end

if rank(obsv(proc_est_)) < size(proc_est_.A, 1)
    disp('Plant is not observable');
else
    disp('Plant is observable');
end

% [~, T] = canon(proc_est_, 'companion');


% proc_est_ = canon(proc_est_, "companion");
% proc_est = ss(proc_est_.A', proc_est_.C', proc_est_.B', proc_est_.D', dt);

% don't need to convert to canonical form since all states are controllable/observable
proc_est = c2d(proc_est_, dt); % forward euler discretization (default)

proc_est

%% tf estimation - 1 pole
Options = tfestOptions;
Options.WeightingFilter = [0 31.4159];
Options.EnforceStability = false;
Options.OutputWeight = [1 0 0;0 1 0;0 0 1];

tf_ = tfest(sys_traindf, [1 1 1;1 1 1;1 1 1], [0 0 0;0 0 0;0 0 0], Options, 'Ts', dt)

%% tf estimation - 2 pole
Options = tfestOptions;
Options.WeightingFilter = [0 31.4159];
Options.EnforceStability = false;
Options.OutputWeight = [1 0 0;0 1 0;0 0 1];

tf_ = tfest(sys_traindf, [2 2 2;2 2 2;2 2 2], [0 0 0;0 0 0;0 0 0], Options, 'Ts', 0.025)

%% convert to state-space
tf_est_ = minreal(ss(tf_), dt);

if rank(ctrb(tf_est_)) < size(tf_est_.A, 1)
    disp('Plant is not controllable');
else
    disp('Plant is controllable');
end

if rank(obsv(tf_est_)) < size(tf_est_.A, 1)
    disp('Plant is not observable');
else
    disp('Plant is observable');
end

% [~, T] = canon(proc_est_, 'companion');


% proc_est_ = canon(proc_est_, "companion");
% proc_est = ss(proc_est_.A', proc_est_.C', proc_est_.B', proc_est_.D', dt);

% don't need to convert to canonical form since all states are controllable/observable
tf_est_.D = zeros(3,3);
tf_est = tf_est_;

tf_est

%% proc estimation
Opt = procestOptions;
Opt.WeightingFilter = [0 31.4159];
Opt.Regularization.Lambda = 10;
proc_ = procest(sys_valdf,{'P1IZ', 'P1IZ', 'P1IZ'; 'P1IZ', 'P1IZ', 'P1IZ'} , Opt);

%% State space model estimation - sim state1/state2: 2 state N4SID
Options = n4sidOptions;
Options.WeightingFilter = [0 31.4159];
Options.Focus = 'simulation';                            
Options.OutputWeight = [1 0 0;0 1 0;0 0 1];              
Options.N4Horizon = [15 15 15];

ss_est = n4sid(sys_traindf, 2, 'Form', 'canonical', Options)

%% State space model estimation - state0: 1 state PEM (time domain)
Options = ssestOptions;
% Options.WeightingFilter = [0 31.4159];
Options.EnforceStability = false;

ss_est = ssest(sys_traind, 1, 'Form', 'canonical', 'DisturbanceModel', 'none', 'Ts', 0.025, Options)

%% State space model estimation - state0: 1 state PEM (freq domain)
Options = ssestOptions;
Options.WeightingFilter = [0 31.4159];
Options.Focus = 'simulation';
Options.OutputWeight = 1;

ss_est = ssest(sys_traindf, 1, 'Form', 'canonical', 'Ts', 0.025, Options)

% make u2 and u3 symmetric
ss_est.B(2) = ( ss_est.B(2) + ss_est.B(3) ) / 2;
ss_est.B(3) = ss_est.B(2);

%% State space model estimation - 2 states
% Options = ssestOptions;
% Options.Focus = 'simulation';
% Options.OutputWeight = [1 0;0 1];
% Options.N4Weight = 'CVA';
% Options.N4Horizon = [15 27 27];

% ss_est = ssest(sys_valdf, 2, 'Form', 'canonical', 'Ts', 0.025, Options)

Options = ssestOptions;
Options.WeightingFilter = [0 31.4159];
Options.Focus = 'simulation';
Options.Regularization.Lambda = 1;
Options.OutputWeight = [1 0;0 1];
Options.SearchOptions.MaxIterations = 50;
Options.N4Horizon = [15 15 15];

ss_est = ssest(sys_valdf, 2, 'Form', 'canonical', 'Ts', 0.025, Options)

%% State space model estimation - 2 states
% Options = ssestOptions;
% % Options.WeightingFilter = [0 31.4159];
% Options.Focus = 'simulation';
% Options.OutputWeight = [1 0;0 1];
% Options.SearchOptions.MaxIterations = 50;
% Options.N4Horizon = [15 15 15];

% ss_est = ssest(sys_valdf, 2, 'Form', 'canonical', 'Ts', 0.025, Options)

Options = ssestOptions;
Options.Focus = 'simulation';
Options.InitialState = 'estimate';
Options.OutputWeight = [1 0;0 1];
Options.SearchOptions.Tolerance = 0;
Options.SearchOptions.MaxIterations = 50;
Options.N4Horizon = [15 15 15];

ss_est = ssest(sys_traindf, 2, 'Form', 'canonical', 'Ts', 0.025, Options)

%% select model
idMdl = ss_est;

%% model validation (prediction error)
figure
compare(idMdl, sys_vald)
print(rootPath + "plots/compare", "-dpng");
rootPath + "plots/compare" + ".png";

%% model validation (residual analysis)
figure
resid(idMdl, sys_vald)
print(rootPath + "plots/resid", "-dpng");
rootPath + "plots/resid" + ".png";

%% augment model
G = idMdl;

if strcmp(idMdlName, 'G0')
    C = [G.C; zeros(2, 1)];
    D = [G.D; zeros(2, 3)];
elseif strcmp(idMdlName, 'G1')
    C = [zeros(1, size(G.C, 2)); G.C];
    D = [zeros(1, 3); G.D];
elseif strcmp(idMdlName, 'G2')
    C = [G.C; zeros(1, 2)];
    D = [G.D; zeros(1, 3)];
end

G = ss(G.A, G.B, C, D, dt);

%% export identified model
save(idMdlName,'G');
