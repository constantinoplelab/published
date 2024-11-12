function [congruent_preferred, congruent_nonpreferred, ...
    incongruent_preferred, incongruent_nonpreferred, info] = ...
    incongruentRewardResponse(SU, S, alignto, zarg)

nAlign = length(alignto);
nc = length(SU);
zscore = @(x) (x - mean(x(:),'omitnan')) / std(x,0,'all','omitnan');

w = arrayfun(@(x) find(round(SU(1).xvec.COFF, 3) == x), [0 0.5]); %using [0 0.5] as wndw - to calculate block significance

[ir, cr] = findIncongruentTrials(S); %find the indices of the first incongruent and congruent trials after a block transition

incongruent_preferred = cell(1, length(alignto));
incongruent_preferred(:) = {nan(length(SU), size(SU(1).hmat.COFF, 2))};
incongruent_nonpreferred = incongruent_preferred;

congruent_preferred = incongruent_preferred;
congruent_nonpreferred = incongruent_preferred;

id = arrayfun(@(x) SU(x).cluster_id, 1:nc)';
sigTable = zeros(length(SU), nAlign);
        
count = zeros(1, nAlign);
for ii = 1:nc
    if zarg
        zSU = arrayfun(@(x) zscore(SU(ii).hmat.(alignto{x})), 1:nAlign,...
            'uniformoutput', false);
    else
        zSU = arrayfun(@(x) SU(ii).hmat.(alignto{x}), 1:nAlign,...
            'uniformoutput', false); %try with raw fr
    end

    for jj = 1:nAlign

        byBlock = arrayfun(@(x) mean(zSU{jj}(S.Block == x, ...
            w(1):w(2)), 2, 'omitnan'), [2, 3], 'UniformOutput', false);

        h = ttest2(byBlock{1}, byBlock{2});

        if h == 1
            sigTable(ii,jj) = 1; %mark if cell has significant block encoding for given alignment

            count(1,jj) = count(1,jj)+1;
            means = arrayfun(@(x) mean(byBlock{x}, 'omitnan'), 1:2);
            if means(1) > means (2)
                pref = 2; %prefers high block; high blocks are marked as 2 in S.Block
            else
                pref = 3; %prefers low block; low blocks are marked as 3 in S.Block
            end
        else
            continue
        end

        if ~isempty(ir.ltom)
            if pref == 2
                incongruent_nonpreferred{jj}(ii, :) = mean(zSU{jj}(ir.ltom(:,1), :), 1, 'omitnan');
            else
                incongruent_preferred{jj}(ii, :) = mean(zSU{jj}(ir.ltom(:,1), :), 1, 'omitnan');
            end
        end

        if ~isempty(ir.htom)
            if pref == 2
                incongruent_preferred{jj}(ii, :) = mean(zSU{jj}(ir.htom(:,1), :), 1, 'omitnan');
            else
                incongruent_nonpreferred{jj}(ii, :) = mean(zSU{jj}(ir.htom(:,1), :), 1, 'omitnan');
            end
        end

        if ~isempty(cr.ltom)
            if pref == 2
                congruent_nonpreferred{jj}(ii, :) = mean(zSU{jj}(cr.ltom(:,1), :), 1, 'omitnan');
            else
                congruent_preferred{jj}(ii, :) = mean(zSU{jj}(cr.ltom(:,1), :), 1, 'omitnan');
            end
        end

        if ~isempty(cr.htom)
            if pref == 2
                congruent_preferred{jj}(ii, :) = mean(zSU{jj}(cr.htom(:,1), :), 1, 'omitnan');
            else
                congruent_nonpreferred{jj}(ii, :) = mean(zSU{jj}(cr.htom(:,1), :), 1, 'omitnan');
            end
        end

    end
end

names = [{'id'} alignto];
sigTable = [id sigTable];
info = array2table(sigTable, 'VariableNames', names);
