function Dat = removeOptOut(Dat)
%function to remove optouts easily from data, keeping indices of the S
%struct and handles consistent

%not exhaustive for masking Dat data. currently only converts what is
%needed for GLM fits and parsing stimuli

%find hits
[~, ~, hits, ~] = parse_choices(Dat.S);
prev_hits = [nan;hits(1:end-1)];

mask = ~isnan(hits);

Dat.hits = hits(mask);
Dat.prev_hits = prev_hits(mask);

%things to change in S from parse_stimuli
if isfield(Dat,'y') %an old form of optimization. eventaully remove
    Dat.y = Dat.y(mask,:);
else
    %new form. update trial numbers
    trial_idx = 1:numel(hits);
    Dat.trial_idx = trial_idx(mask);    
end
Dat.S.pd{1}.right_prob = Dat.S.pd{1}.right_prob(mask);
Dat.S.pd{1}.left_prob = Dat.S.pd{1}.left_prob(mask);
Dat.S.pd{1}.went_right = Dat.S.pd{1}.went_right(mask);
Dat.S.pd{1}.went_left = Dat.S.pd{1}.went_left(mask);

%flatten cells to long array for processing via binspikes
%lflash_vec = [Dat.handles.lflashes{:}];
%rflash_vec = [Dat.handles.lflashes{:}];
%lbups_vec = [Dat.handles.lbups{:}];
%rbups_vec = [Dat.handles.rbups{:}];

%handle Dat.handles
Dat.handles.start = Dat.handles.start(mask);
Dat.handles.end = Dat.handles.end(mask);
Dat.handles.lflashes = Dat.handles.lflashes(mask);
Dat.handles.rflashes = Dat.handles.rflashes(mask);
Dat.handles.lbups = Dat.handles.lbups(mask);
Dat.handles.rbups = Dat.handles.rbups(mask);
Dat.handles.went_right = Dat.handles.went_right(mask);
Dat.handles.choice = Dat.handles.choice(mask);
Dat.handles.leavecpoke = Dat.handles.leavecpoke(mask);

%S struct for parse choices
%find values or probability and water volume, case of rat choosing right port
Dat.S.pd{1}.hits = Dat.S.pd{1}.hits(mask);
Dat.S.pd{1}.this_right_volume = Dat.S.pd{1}.this_right_volume(mask);
Dat.S.pd{1}.this_left_volume = Dat.S.pd{1}.this_left_volume(mask);
Dat.S.pd{1}.right_hit = Dat.S.pd{1}.right_hit(mask);
Dat.S.pd{1}.left_hit = Dat.S.pd{1}.left_hit(mask);

