% postprocessing for all ofc data (stratified and non-stratified)


%location of analysis codebase
codepath = '~/projects/constantinoplelab/Analysis/';

% add code paths
addpath(genpath(codepath))
addpath(genpath(codepath + "david"))
addpath(genpath(codepath + "maggie"))

%for saving pds correctly so illustrator can modify them
set(0, 'DefaultFigureRenderer', 'painters');

% DS projection neuron data location
datadir = "/Users/dhocker/projects/dynamics/data/maggie/";
savedir = "/Users/dhocker/projects/dynamics/results/maggie/";

fname = datadir + "DS_projection_neurons_non-Stimulated.mat";
b = load(fname);


%%
% parse the DS projections to get beliefs
usetorch = false

decoderlist = {'psth2','none'};
%decoderlist = {'psth2'};
epochlist = {'reward','coff','son'} % son didn't get done for some reason
%epochlist = {'son'} % son didn't get done for some reason
usemlblist = {true, false};
%stratifiedlist = {true, false}
stratifiedlist = {false}

for is_stratified = stratifiedlist
    is_stratified = is_stratified{1};
    for decodertype = decoderlist
         decodertype = decodertype{1};
        for epoch = epochlist
            epoch = epoch{1};
            for usemlb = usemlblist
                usemlb = usemlb{1};
    
                disp(epoch)
                disp(decodertype)
                disp(usemlb)
    
                output = dSTR_decode_parsedata(usetorch, epoch,decodertype);
                
                if usemlb
                    if is_stratified
                        dfun = @(x) strcat(savedir,'allofc_stratified/','mlb_allofc_',num2str(x),'_',epoch,'_',decodertype,'_stratified_',num2str(x),'.mat');
                    else
                        dfun = @(x) strcat(savedir,'allofc/','mlb_allofc_',num2str(x),'_',epoch,'_',decodertype,'_',num2str(x),'.mat');
                    end
                else
                    if is_stratified
                        dfun = @(x) strcat(savedir,'allofc_stratified/','mlb2_allofc_',num2str(x),'_',epoch,'_',decodertype,'_stratified_',num2str(x),'.mat');
                    else
                        dfun = @(x) strcat(savedir,'allofc/','mlb2_allofc_',num2str(x),'_',epoch,'_',decodertype,'_',num2str(x),'.mat');
                    end
                end
                
                
                % load the data and get preds and true for each fold
                %nshuff = 20;
                nshuff = 100;
                nfold = 20;
                xvec_ofc = -4:0.2:4;
                ntx_ofc = numel(xvec_ofc);
                ntrials_test = 300;
                
                d = struct();
                d.preds_CV = nan(nshuff, nfold, ntx_ofc, ntrials_test);
                d.true_CV = nan(nshuff, nfold, ntx_ofc, ntrials_test);
                
                for j = 1:nshuff
                    %disp(j)
                    f = load(dfun(j));
                    d.preds_CV(j,:,:,:) = f.preds_CV;
                    d.true_CV(j,:,:,:) = f.true_CV;
                end
                
                % calculate the accuracy
                
                
                
                preds_all_ofc = zeros(nshuff,nfold,ntx_ofc);
                for n = 1:nshuff
                    for m = 1:nfold
                        for j = 1:ntx_ofc
                            preds_all_ofc(n,m,j) = sum(d.preds_CV(n,m,j,:)==d.true_CV(n,m,j,:))/numel(d.true_CV(n,m,j,:));
                        end
                    end
                end
                
                accuracy_pred_means_ofc = squeeze(mean(preds_all_ofc,2,'omitnan'));
                accuracy_ofc_mean = mean(accuracy_pred_means_ofc,1,'omitnan');
                accuracy_ofc_median = median(accuracy_pred_means_ofc,1,'omitnan');
                accuracy_ofc_sem = std(accuracy_pred_means_ofc,[],1,'omitnan')/sqrt(nshuff); %number of shuffles
                
                
                % save
                
                if usemlb
                    if is_stratified
                        savename =strcat(savedir,'allofc_stratified/','postprocess_mlb_allofc_',epoch,'_',decodertype,'_stratified.mat');
                    else
                        savename = strcat(savedir,'allofc/','postprocess_mlb_allofc_',epoch,'_',decodertype,'.mat');
                    end
                else
                    if is_stratified
                        savename = strcat(savedir,'allofc_stratified/','postprocess_mlb2_allofc_',epoch,'_',decodertype,'_stratified.mat');
                    else
                        savename = strcat(savedir,'allofc/','postprocess_mlb2_allofc_',epoch,'_',decodertype,'.mat');
                    end
                end
                
                disp(savename)
                save(savename, 'accuracy_ofc_mean', 'accuracy_ofc_sem','accuracy_ofc_median','preds_all_ofc')
            end
        end
    end
end


