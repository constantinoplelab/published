function plot_fig1f(datadir)
% Plot event-aligned head speed and head angle relative to the central nose
% port for an example session.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.


load(fullfile(datadir, 'data-published', ...
    'J063_AChDA_20230201_HJJ_DLC.mat'), 'AlignmentStruct')
plotDLCaligned(AlignmentStruct);


end