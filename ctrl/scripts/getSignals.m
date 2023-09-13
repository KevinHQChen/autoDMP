function signals = getSignals(sim, out)
    warning('off', 'SimulationData:Objects:InvalidAccessToDatasetElement');
    signalNames = ["t"];
    signalValues = [out.tout];
    signals = containers.Map(signalNames, signalValues);

    signalHandles = find_system(sim,'FindAll','on','type','line');
    for i = 1:size(signalHandles, 1)
        % find names of all signals in simulation
        newSigName = string(get_param(signalHandles(i), 'Name'));
        % discard unnamed, repeated, and unlogged signals
        if ~isempty(newSigName) & newSigName ~= "" & ~ismember(newSigName, signalNames) & ~isempty(get(out.logsout, newSigName))
            signalNames = [signalNames newSigName];
            signals(newSigName) = out.logsout.get(newSigName).Values.Data;
        end
    end
end

% function signals = getSignals(sim, t, out)
%     warning('off', 'SimulationData:Objects:InvalidAccessToDatasetElement');
%     signalNames = ["t"];
%     signals = containers.Map(signalNames, t);

%     signalHandles = find_system(sim,'FindAll','on','type','line');
%     for i = 1:size(signalHandles, 1)
%         % find names of all signals in simulation
%         newSigName = string(get_param(signalHandles(i), 'Name'));
%         % discard unnamed, repeated, and unlogged signals
%         if ~isempty(newSigName) & newSigName ~= "" & ~ismember(newSigName, signalNames) & ~isempty(get(out.logsout, newSigName))
%             signalNames = [signalNames newSigName];
%             signals(newSigName) = out.logsout.get(newSigName).Values.Data;
%         end
%     end
% end
