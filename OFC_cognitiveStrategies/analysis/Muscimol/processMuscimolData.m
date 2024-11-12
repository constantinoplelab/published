function [control, muscimol] = processMuscimolData(ratList1, ratList2, ...
    muscimolDataPath, twin, binSize, smoothfactor)

nrats = length(ratList1);

nback = 10;

%preallocate matrices
betasC = nan(nrats, nback+1);
betasM = betasC;

wtPrevRewC = nan(nrats, 2);
wtPrevRewM = wtPrevRewC;


for rr = 1:nrats

    S1 = load(strcat(muscimolDataPath, ratList1{rr}, '.mat'));
    S2 = load(strcat(muscimolDataPath, ratList2{rr}, '.mat'));
    
    %wait time curves
    [hiC(rr), loC(rr), mixC(rr)] = wtcurves(S1.A);
    [hiE(rr), loE(rr), mixE(rr)] = wtcurves(S2.A);
    
    %regress wait time vs reward
    [control.regress(rr,:), ~] = regress_wt_vs_rew(S1.A, nback);
    [muscimol.regress(rr,:), ~] = regress_wt_vs_rew(S2.A, nback);

    %conditional wait time for 20ul
    [control.prevRew(rr,:), ~] = waittime_20ul_by_previous_vol(S1.A);
    [muscimol.prevRew(rr,:), ~] = waittime_20ul_by_previous_vol(S2.A);
    
    %wait time dynamics
    [~, ~, control.mtol(rr,:), control.mtoh(rr,:), control.ltom(rr,:), control.htom(rr,:), ~] =...
        block_dynamics_wt_binTrials(S1.A, twin, binSize, smoothfactor);

    [~, ~, muscimol.mtol(rr,:), muscimol.mtoh(rr,:), muscimol.ltom(rr,:), muscimol.htom(rr,:), ~] =...
        block_dynamics_wt_binTrials(S2.A, twin, binSize, smoothfactor);

end

%wait time curves
control.hi = cell2mat(arrayfun(@(x) hiC(x).wt, 1:nrats, 'uniformoutput', false)');
control.lo = cell2mat(arrayfun(@(x) loC(x).wt, 1:nrats, 'uniformoutput', false)');
control.mix = cell2mat(arrayfun(@(x) mixC(x).wt, 1:nrats, 'uniformoutput', false)'); 

muscimol.hi = cell2mat(arrayfun(@(x) hiE(x).wt, 1:nrats, 'uniformoutput', false)');
muscimol.lo = cell2mat(arrayfun(@(x) loE(x).wt, 1:nrats, 'uniformoutput', false)');
muscimol.mix =  cell2mat(arrayfun(@(x) mixE(x).wt, 1:nrats, 'uniformoutput', false)');
