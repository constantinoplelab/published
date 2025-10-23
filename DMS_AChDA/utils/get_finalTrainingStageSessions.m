function S_final = get_finalTrainingStageSessions(S)
% get training stage 9 days in S struct and remove everything else
usethese = cellfun(@(s) isstruct(s), S.pd);

S_final = S;
S_final.pd = S.pd(usethese);
S_final.peh = S.peh(usethese);

stage = cell2mat(cellfun(@(sess) sess.TrainingStage(1), S_final.pd',...
    UniformOutput=false));

S_final.pd(stage<9) = [];
S_final.peh(stage<9) = [];

end