function [d_fr, d_iti, idx] = get_dFRdITI(neuron, S, A)

% look twin sec around the event
twin = [0 0.3];
T = neuron.xvec.(A{1});
[~,t1] = min(abs(T-twin(1)));
[~,t2] = min(abs(T-twin(2)));

% filter out outliers in ITI
cutoff = [prctile(S.iti,2), prctile(S.iti,98)];
fprintf('ITI cutoff (s): %.2f, %.2f\n', cutoff)
ITI = S.iti;
ITI(ITI<min(cutoff) | ITI>max(cutoff)) = nan;

for a=1:length(A)

    % compares between closest trial of the same type and ignore other
    % trials in between
    if strcmpi(A{a}, 'COFF') % all trials
        idx{a} = [1:length(S.Block)]';
    elseif strcmpi(A{a}, 'SON') % non-vios
        idx{a} = find(~S.vios);
    elseif sum(strcmpi(A{a}, {'SOFF', 'Rew'}))>0 % hit trials
        idx{a} = find(S.hits);
    elseif strcmpi(A{a}, 'Opt') % opt-outs
        idx{a} = find(S.optout);
    end

    data = neuron.hmat.(A{a})(idx{a}, t1:t2);
    avgfr = mean(data, 2, 'omitnan');
    d_fr{a} = diff(avgfr);
    d_iti{a} = -log2(ITI(idx{a}(2:end))) + log2(ITI(idx{a}(2:end)-1));

end


end


