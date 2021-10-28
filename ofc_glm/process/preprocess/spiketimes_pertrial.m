function spikes_pertrial = spiketimes_pertrial(spiketimes, handles, wndw)
%SPIKETIMES_PERTRIAL parses spike timings into different trials
%
%output: 
%spikes_pertrial: 

nt = length(handles.start);
spikes_pertrial = cell(nt,1);

%find trial spikes by searchign range between trial start and end/violation
for j = 1:nt
    %if rat participates in trial use handles.end
    if isfield(handles, 'end')
        these = find(spiketimes>=handles.start(j)-wndw & spiketimes<=handles.end(j)+wndw);
        spikes_pertrial{j} = spiketimes(these)-handles.start(j);
    %if not, use the time of leaving center poke
    else
        these = find(spiketimes>=handles.start(j)-wndw & spiketimes<=handles.leavecpoke(j)+wndw);
        spikes_pertrial{j} = spikestimes(these)-handles.start(j);
    end
end