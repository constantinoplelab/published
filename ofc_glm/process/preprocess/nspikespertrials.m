function n = nspikespertrials(spiketimes, handles, wndw)
%NSPIKESPERTRIALS Counts number or spikes on each trial in a session
%   used to determine if neuron is active enough to use for analysis
%   returns (n x 1) array with number of spikes, where n = number trials

n = nan(length(handles.start), 1);

%find trial spikes by searchign range between trial start and end/violation
for j = 1:length(handles.start)
    %if rat participates in trial use handles.end
    if isfield(handles, 'end')
        these = find(spiketimes>=handles.start(j)-wndw & spiketimes<=handles.end(j)+wndw);
    %if not, use the time of leaving center poke
    else
        these = find(spiketimes>=handles.start(j)-wndw & spiketimes<=handles.leavecpoke(j)+wndw);
    end
    n(j) = numel(these);
end