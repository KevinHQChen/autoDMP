clc, close all, clear all
rootPath = "~/autoDMP/ctrl/scripts/";
dataset = 'train'; % [train | val]
plantConstants;

%% supervisor dependencies
nullEv = struct('r', [0; 0; 0; 0; 0; 0], 'preT', 0, 'moveT', 0, 'postT', 0);
busInfo = Simulink.Bus.createObject(nullEv);
event_bus = evalin('base', busInfo.busName);
evalin('base', ['clear ' busInfo.busName]);

% mpc stuff
beta = 0.13534; % maximize robustness in closed-loop performance
ywt0 = zeros(1, 2*no);

%% generate excitation signal
% PRBS
range = [-1, 1];
clockPeriod = 4;
band = [0 1/clockPeriod];
u_id = idinput([clockPeriod*2^10-3-1-1, 3], 'prbs', band, range);
% swap inputs to minimize correlation btwn train/val datasets
if strcmp(dataset, 'train')
    u_id = [u_id(:,1), u_id(:,2), u_id(:,3)];
else
    u_id = [u_id(:,3), u_id(:,1), u_id(:,2)];
end
u_id = u_id + repmat(u_0', size(u_id, 1), 1);  % laplace pressure compensation
u_id = iddata([], u_id, dt);

u = u_id.InputData;
t = u_id.SamplingInstants;

figure;
plot(u_id);
grid
figure;
plot(fft(u_id));
grid
figure;
[pxx, f] = pwelch(u_id.InputData, [], [], [], samplingRate);
plot(f, 10*log10(pxx));
grid
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density of PRBS Signal');

%% run simulation and analyze results
sim = 'T_junction_mdl';
open(sim);
commandwindow;

% input data setup
% WARNING for multidimensional inputs we have to manually set inport dimensions
set_param(sim, 'LoadExternalInput', 'on',...
               'ExternalInput', '[t, u]');

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
% stepSize = 1e-3;
% set_param(sim, 'StopTime', num2str(t(end)),...
%                'SolverType', 'Fixed-step',...
%                'Solver', 'FixedStepAuto',...
%                'FixedStep', num2str(stepSize));
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
% plotSignals(signals);

%% store results in file
ySignal = squeeze(signals('y'))';

writetable(...
    array2table([signals('u'), ySignal(:, 1:3)],...
                'VariableNames',...
                {'u0', 'u1', 'u2', 'y0', 'y1', 'y2'}),...
    rootPath + "simSysID/" + dataset + "/ctrlDataQueue.txt");
disp('input-output data stored in ' + rootPath + "simSysID/" + dataset + "/ctrlDataQueue.txt");
