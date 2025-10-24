function [avgFR, RPE_firstTrial, firstTrialIdx] = ...
    get_avgFR_persistence(neuron, beh, A, twin)
% get average firing rate at A on a sequence of 3 trials where the RPE
% on the 2nd trial is small 

if nargin<4
    twin = [0 0.3];
end

% get an RPE sequence    
% filter out outliers in ITI and compute -delta ITI
cutoff = [prctile(beh.iti,2), prctile(beh.iti,98)];
iti = beh.iti;
iti(iti<min(cutoff) | iti>max(cutoff)) = nan;
dITI = -diff(log2(iti));

% get |dITI|<0.2 trials
thresh = 0.2;
smallRPEtrials = find(abs(dITI)<thresh);
% If dITI(N) = -(ITI(N+1)-ITI(N)) ~ 0, means RPE on Nth trial ~ 0

spikes = neuron.hmat.(A);
FRs = cell(length(smallRPEtrials),1);
RPE_firstTrial = nan(length(smallRPEtrials),1); % RPE on trial N-1

for tr=1:length(smallRPEtrials)
    N = smallRPEtrials(tr);
    % get spikes
    FRs{tr} = spikes(N-1:N+1,:);
    % get RPE on trial N-1
    RPE_firstTrial(tr) = dITI(N-1);
end

validTrials = find(cellfun(@(x) ~isempty(x), FRs));

if ~isempty(validTrials)
    T = neuron.xvec.COFF;
    [~,t1] = min(abs(T-twin(1)));
    [~,t2] = min(abs(T-twin(2)));

    avgFR = nan(length(validTrials),3);
    for tr=1:length(validTrials)
        data = FRs{validTrials(tr)};
        avgFR(tr,:) = mean(data(:,t1:t2), 2, 'omitnan')';
    end
    RPE_firstTrial = RPE_firstTrial(validTrials);
end

% index of the first trial of the sequence (N-1) 
firstTrialIdx = smallRPEtrials(validTrials)-1;

end