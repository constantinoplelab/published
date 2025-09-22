function [DA_offercue, DA_sideon, DA_rewardcue, DA_reward, DA_optout] = ...
    get_avgFl_allevents(datadir, ratlist, sensor, region)

basedir = fullfile(datadir, 'data-published/PhotometryDataAligned');
pdatadir = fullfile(basedir, strcat(sensor, '_', region));

rewards = [5, 10, 20, 40, 80];
DA_offercue = cell(1, length(rewards));

DA_sideon = struct;
DA_sideon.contra = [];
DA_sideon.ipsi = [];
DA_optout = struct;
DA_optout.contra = [];
DA_optout.ipsi = [];

nbin=4;
DA_rewardcue = cell(1, nbin);
DA_reward = cell(1, nbin);

for rr=1:length(ratlist)
    rat = ratlist{rr};

    load(fullfile(pdatadir, strcat(rat, '_avgFbyVol_bc.mat')), 'data');
    for rew=1:length(rewards)
        DA_offercue{rew}(rr,:) = data{rew,2}; % coff=2
    end

    load(fullfile(pdatadir, strcat(rat, '_avgFbySide_bc.mat')), ...
        'dataContra', 'dataIpsi');
    DA_sideon.contra(rr,:) = dataContra{3}; % sideon=3
    DA_sideon.ipsi(rr,:) = dataIpsi{3};
    DA_optout.contra(rr,:) = dataContra{6}; % optout=6
    DA_optout.ipsi(rr,:) = dataIpsi{6};

    load(fullfile(pdatadir, strcat(rat, '_avgFbyDelay_bc.mat')), 'data');
    for bb=1:nbin
        DA_rewardcue{bb}(rr,:) = data{bb,4}; % soff=4
        DA_reward{bb}(rr,:) = data{bb,5}; % reward=5
    end

end


end



