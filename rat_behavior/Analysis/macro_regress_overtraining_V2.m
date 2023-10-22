function [B1, B2, B1iti, B2iti,...
    dwt1, dwt2, diti1, diti2,...
    pcatch1, pcatch2, pvio_iti1, pvio_iti2] =...
    macro_regress_overtraining_V2(nsess, ratList, Sdir)

%this version compares early sessions in training stage 9 from the S
%struct, instead of the A struct, which has already been curated. cmc
%11/16/22.
% updated to incorporate into Figure 3 code. avm 27 Jan 2023
grp = 2;

[B1, B2, B1iti, B2iti] = deal(cell(1, 3));

for k = 1:3
    B1{k} = nan(length(ratList), ceil(nsess/grp));
    B2{k} = nan(length(ratList), ceil(nsess/grp));
    B1iti{k} = nan(length(ratList), ceil(nsess/grp));
    B2iti{k} = nan(length(ratList), ceil(nsess/grp));
end

dwt1 = nan(length(ratList), ceil(nsess/grp)); dwt2 = dwt1;
diti1 = dwt1; diti2 = dwt1;
pcatch1 = dwt1; pcatch2 = dwt1;
pvio_iti1 = dwt1;
pvio_iti2 = dwt2;


for rr = 1:length(ratList)
    disp(rr);

    % Load Data
    a = load([Sdir 'ratTrial_' ratList{rr} '.mat']);
    A = a.A;
    
    % fname = strcat(['SStruct_', ratList{rr}, '.mat']);
    % s = load(strcat([Sdir, fname])); S = s.S;
    % A = parse_data_from_mysql(S, [], 'all');
    
    A.wait_time(A.wait_time>A.wait_thresh) = nan;
    A.ITI(A.ITI>prctile(A.ITI,99)) = nan;
    ctr = 1;

    firstn = 1:nsess;
    lastn = length(A.ntrials)-nsess+1:length(A.ntrials);

    tctr1 = 1;
    tctr2 = sum(A.ntrials(1:length(A.ntrials)-nsess))+1;
    for k = 1:grp:nsess

        these = firstn(k):firstn(k)+grp-1;
        if sum(these>length(A.ntrials))>0 %if indices surpass data, nan it.
            beta1 = nan(3,1);
            beta1iti = nan(3,1);
            dwt1(rr,ctr) = nan;
            diti1(rr,ctr) = nan;
            pcatch1(rr,ctr) = nan;
            pvio_iti1(rr,ctr) = nan;
        else
            idx = tctr1:tctr1+sum(A.ntrials(these))-1;
            [A1] = parseSstruct(A, idx, 0);

            tctr1 = tctr1+sum(A.ntrials(these));
            [beta1, beta1iti, dwt1(rr,ctr), diti1(rr,ctr)] =...
                regress_trials_blocksv2(A1, 1);
            pcatch1(rr,ctr) =...
                mean(A1.optout, 'omitnan')/mean(A1.prob_catch, 'omitnan');

            [A1] = parseSstruct(A, idx, 1);
            [~, ~, ~, ~, pvio_iti1(rr,ctr)] =...
                regress_trials_blocksv2(A1, 1);
        end

        these = lastn(k):lastn(k)+grp-1;
        if sum(these<=0)>0 %if there are negative indices, nan it.
            beta2 = nan(3,1);
            beta2iti = nan(3,1);
            dwt2(rr,ctr) = nan;
            diti2(rr,ctr) = nan;
            pcatch2(rr,ctr) = nan;
            pvio_iti2(rr,ctr) = nan;
        else

            idx = tctr2:tctr2+sum(A.ntrials(these))-1;
            [A2] = parseSstruct(A, idx, 0);

            tctr2 = tctr2+sum(A.ntrials(these));
            [beta2, beta2iti, dwt2(rr,ctr), diti2(rr,ctr)] =...
                regress_trials_blocksv2(A2, 1);
            pcatch2(rr,ctr) =...
                mean(A2.optout, 'omitnan')/mean(A2.prob_catch, 'omitnan');

            [A2] = parseSstruct(A, idx, 1);
            [~, ~, ~, ~, pvio_iti2(rr,ctr)] =...
                regress_trials_blocksv2(A2, 1);
        end
        for jj = 1:3
            B1{jj}(rr,ctr) = beta1(jj);
            B2{jj}(rr,ctr) = beta2(jj);

            B1iti{jj}(rr,ctr) = beta1iti(jj);
            B2iti{jj}(rr,ctr) = beta2iti(jj);
        end
        ctr = ctr+1;
    end

    dwt1(rr,:) = smooth(dwt1(rr,:));
    dwt2(rr,:) = smooth(dwt2(rr,:));

    diti1(rr,:) = smooth(diti1(rr,:));
    diti2(rr,:) = smooth(diti2(rr,:));

    B1{jj}(rr,:) = smooth(B1{jj}(rr,:));
    B2{jj}(rr,:) = smooth(B2{jj}(rr,:));

    B1iti{jj}(rr,:) = smooth(B1iti{jj}(rr,:));
    B2iti{jj}(rr,:) = smooth(B2iti{jj}(rr,:));

    pcatch1(rr,:) = smooth(pcatch1(rr,:));
    pcatch2(rr,:) = smooth(pcatch2(rr,:));

    pvio_iti1(rr,:) = smooth(pvio_iti1(rr,:));
    pvio_iti2(rr,:) = smooth(pvio_iti2(rr,:));

end

dwt1 = rmvoutlier(dwt1);
dwt2 = rmvoutlier(dwt2);
diti1 = rmvoutlier(diti1);
diti2 = rmvoutlier(diti2);
pcatch1 = rmvoutlier(pcatch1);
pcatch2 = rmvoutlier(pcatch2);
pvio_iti1 = rmvoutlier(pvio_iti1);
pvio_iti2 = rmvoutlier(pvio_iti2);
end

function [A1] = parseSstruct(A, ix, usepvio)

A1.block = A.block(ix);
A1.wait_time = A.wait_time(ix);
A1.reward = A.reward(ix);
A1.vios = A.vios(ix);
A1.optout = A.optout(ix);
A1.hits = A.hits(ix);
A1.ITI = A.ITI(ix);
A1.prob_catch = A.prob_catch(ix);

if usepvio==0
    bad = find(A1.vios==1);
    A1.block(bad) = [];
    A1.wait_time(bad) = [];
    A1.reward(bad) = [];
    A1.vios(bad) = [];
    A1.optout(bad) = [];
    A1.hits(bad) = [];
    A1.ITI(bad) = [];
    A1.prob_catch(bad) = [];
elseif usepvio==1
    bad = find(A1.vios~=1);
    bad = bad+1; bad(bad>length(A1.vios)) = [];
    A1.block(bad) = [];
    A1.wait_time(bad) = [];
    A1.reward(bad) = [];
    A1.vios(bad) = [];
    A1.optout(bad) = [];
    A1.hits(bad) = [];
    A1.ITI(bad) = [];
    A1.prob_catch(bad) = [];
end
end

function [x] = rmvoutlier(x)
x(isoutlier(x, 'movmedian', 5)) = nan;
end
