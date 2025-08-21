function ProcessData_ExtendedData9(datadir, savedir, codedir)
%ProcessData_ExtendedData9 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData9.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_ExtendedData9.mat' from Zenodo was saved
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved
%   savedir - Local directory where you would like the outputs to be saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it isn't already
    addpath(genpath(codedir))
end

%% Load raw data
load([datadir, 'RawData_ExtendedData9.mat'], 'PvsD_DOWN_NAccbackground', ...
    'EvsD_DOWN_NAccbackground', 'EvsD_UP_NAccbackground', ...
    'dopaminegenes_lists', 'PvsD_foldchange', 'EvsD_foldchange')

%% Process data
%--------------------------------------------------------------------------
%ED9a. Proestrus vs diestrus: downregulated genes compared to
% all identified proteins in NAcc as background, using STRING v12
%--------------------------------------------------------------------------
[~,sortedPvsDDOWNidx] = sort(PvsD_DOWN_NAccbackground.('strength'));
sortedPvsD_DOWN_NAccbackground = PvsD_DOWN_NAccbackground(sortedPvsDDOWNidx,:);

%--------------------------------------------------------------------------
%ED9b. Estrus vs. diestrus: downregulated proteins compared to
% all identified proteins in NAcc as background, using STRING v12
%--------------------------------------------------------------------------
[~,sortedEvsDDOWNidx] = sort(EvsD_DOWN_NAccbackground.('strength'));
sortedEvsD_DOWN_NAccbackground = EvsD_DOWN_NAccbackground(sortedEvsDDOWNidx,:);

%--------------------------------------------------------------------------
%ED9c. Estrus vs diestrus: upregulated proteins compared to
% all identified proteins in NAcc as background, using STRING v12
%--------------------------------------------------------------------------
[~,sortedEvsDUPidx] = sort(EvsD_UP_NAccbackground.('strength'));
sortedEvsD_UP_NAccbackground = EvsD_UP_NAccbackground(sortedEvsDUPidx,:);

%--------------------------------------------------------------------------
%ED9d. Heatmap of dopamine-related terms, differential expression in
%our samples. Highlighting dopamine reuptake as term to change.
%--------------------------------------------------------------------------
%import dopamine-related proteins from https://geneontology.org/ terms related to dopamine for rattus norvegicus
dopaminegenes = unique(lower(dopaminegenes_lists.GENES));
dopamine_uptake_genes = {'slc6a3','slc6a4','smpd3'};

PvsD_FC = PvsD_foldchange.Log2_Low_Medium_;
PvsD_coverage = PvsD_foldchange.SequenceCoverage___;
PvsD_proteins = lower(PvsD_foldchange.Gene);
PvsD_proteins = PvsD_proteins(strcmp(PvsD_foldchange.Student_sT_testSignificantLow_Medium, '+')); %only keep significant genes

EvsD_FC = EvsD_foldchange.Log2_Low_High_;
EvsD_coverage = EvsD_foldchange.SequenceCoverage___;
EvsD_proteins = lower(EvsD_foldchange.Gene);

PvsD_foldchange_DA = NaN(length(dopaminegenes), 1);
EvsD_foldchange_DA = NaN(length(dopaminegenes), 1);
for g = 1:length(dopaminegenes)
    gene = dopaminegenes{g};
    if sum(strcmp(PvsD_proteins, gene))==1
        PvsD_foldchange_DA(g,1) = PvsD_FC(strcmp(PvsD_proteins, gene));
    elseif sum(strcmp(PvsD_proteins, gene))>1
        multiplegenes = find(strcmp(PvsD_proteins, gene));
        [~, cvrg_idx] = max(PvsD_coverage(multiplegenes));
        PvsD_foldchange_DA(g,1) = PvsD_FC(multiplegenes(cvrg_idx));
    elseif sum(strcmp(PvsD_proteins, gene))==0
        PvsD_foldchange_DA(g,1) = NaN;
    end
    if sum(strcmp(EvsD_proteins, gene))==1
        EvsD_foldchange_DA(g,1) = EvsD_FC(strcmp(EvsD_proteins, gene));
    elseif sum(strcmp(EvsD_proteins, gene))>1
        multiplegenes = find(strcmp(EvsD_proteins, gene));
        [~, cvrg_idx] = max(EvsD_coverage(multiplegenes));
        EvsD_foldchange_DA(g,1) = EvsD_FC(multiplegenes(cvrg_idx));
    elseif sum(strcmp(EvsD_proteins, gene))==0
        EvsD_foldchange_DA(g,1) = NaN;
    end
end
% remove NaNs (same in proestrus and estrus)
PvsD_foldchange_DA_iaN = PvsD_foldchange_DA(~isnan(PvsD_foldchange_DA));
EvsD_foldchange_DA_iaN = EvsD_foldchange_DA(~isnan(EvsD_foldchange_DA));
dopaminegenes_iAN = dopaminegenes(~isnan(PvsD_foldchange_DA));

%sort by fold change
[~, sortidx] = sort(PvsD_foldchange_DA_iaN, 'descend');
sortedPvsD = -1*PvsD_foldchange_DA_iaN(sortidx); %positive=greater in proestrus
sortedDAgenes = dopaminegenes_iAN(sortidx);

[~, sortidx] = sort(EvsD_foldchange_DA_iaN, 'descend');
sortedEvsD = -1*EvsD_foldchange_DA_iaN(sortidx); %positive=greater in proestrus

%% Save processed data 
save([savedir 'ProcessData_ExtendedData9'], ...
    'sortedPvsD_DOWN_NAccbackground', 'sortedEvsD_DOWN_NAccbackground', ...
    'sortedEvsD_UP_NAccbackground', 'dopamine_uptake_genes', ...
    'sortedPvsD', 'sortedDAgenes', 'sortedEvsD', '-v7.3');

end