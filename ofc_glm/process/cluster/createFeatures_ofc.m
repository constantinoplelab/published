function [fet_proc,fet,fet_sem,fet_legend,fet_inds,xvec] =...
        createFeatures_ofc(ops)
    %generates features from OFC data. 
    %The running history of the many
    %types investigated, but pulled and made on the fly
    % includes preprocessing steps such as denoising via PCA
    
    
 
%% load data and sizing

fname = strcat(ops.savedir,ops.fname_base);

f = load(fname);
s = f.sourceops; %the vector information

nbase_EV = numel(f.EV_bins)-1;
nbase_rew = numel(f.EV_bins);



%% which class of features

switch ops.featureClass
    case 'singleConditional'
        switch ops.featureType
                
            case 1 %choice and outcome, expected value L/R presented offer, rewarded volume
                fetCO = cat(2,f.COfet{1:6});
                fetCO_sem = cat(2,f.COfet_sem{1:6});
                fetCO_legend = f.COfet_legend(1:6);
                
                fetPO = cat(2,f.POfet{5:6});
                fetPO_sem = cat(2,f.POfet_sem{5:6});
                fetPO_legend = f.POfet_legend(5:6);
                
                fetreward = f.Rewfet{5};
                fetreward_sem = f.Rewfet_sem{5};    
                fetreward_legend = f.Rewfet_legend(5);
                
                fet = cat(2,fetCO,fetPO,fetreward);
                fet_sem = cat(2,fetCO_sem,fetPO_sem,fetreward_sem);
                fet_legend = [fetCO_legend,fetPO_legend,fetreward_legend];
                fet_inds = [1:7,7+nbase_EV,7+2*nbase_EV,7+2*nbase_EV + nbase_rew];
                
                xvec_CO = {1:6};
                xvec_PO = {s.EV_bins(1:end-1),s.EV_bins(1:end-1)};
                xvec_Rew = {s.EV_bins};
                xvec = {xvec_CO,xvec_PO,xvec_Rew};
             
        end

        
    case 'PSTH'
        switch ops.featureType %aligned to different events
            case 1
                fet = f.psth_fet;
                fet_sem = f.psth_fet_sem;
                fet_legend = {'psth'};
                xvec = {f.tmesh};
                fet_inds = 1;
                
            case 2
                fet = f.psth_fet_leavecpoke;
                fet_sem = f.psth_fet_leavecpoke_sem;
                fet_legend = {'psth'};
                xvec = {f.tmesh_leavecpoke};
                fet_inds = 1;
        end
        
   
end


%%  preprocess 

fet_proc = fet;
nfet = size(fet,2);
nsamp = size(fet,3);

%scrub nans
fet_proc(isnan(fet_proc))=0;
fet(isnan(fet))=0;


%find one PC set based on the full dataset, then use those PC on the rest
%of the samples

fetj = fet_proc(:,:,1);
cutoffval = ops.preprocCutoff;
[coeff_old,score,latent] = pca(fetj);
cutoffvec = find(cumsum(latent)/sum(latent) > cutoffval);
cutoff = cutoffvec(1);

%decide on new sizing
fet_proc_old = fet; %used to projectig old feature space onto common pc
fet_proc = fet(:,1:cutoff,:);
fet_proc(:,:,1) = score(:,1:cutoff);


for j = 2:nsamp
    
    fetj = fet_proc(:,:,j);  
        
    switch ops.preprocType

        case 'svd' %TODO: this is just a desnoiing form I think
            cutoffval = ops.preprocCutoff;              
            [U,S,V] = svd(fetj);
            cutoffvec = find(cumsum(diag(S))/sum(diag(S)) > cutoffval); 
            cutoff = cutoffvec(1);
            fet_proc(:,:,j) = U*S(:,1:cutoff)*V(:,1:cutoff)';       

        case 'pca' %projected onto principal components capturing >95% variance
            fetj = fet_proc_old(:,:,j); 
            score = fetj*coeff_old;
            fet_proc(:,:,j) = score(:,1:cutoff);
            
        case 'pca_nocutoff' %projected onto principal components capturing >95% variance
            cutoffval = ops.preprocCutoff;
            [~,score,latent] = pca(fetj);
            cutoffvec = find(cumsum(latent)/sum(latent) > cutoffval);
            cutoff = cutoffvec(1);
            fet_proc(:,1:cutoff,j) = score(:,1:cutoff);
            fet_proc(:,cutoff+1:nfet,j) = 0; %if feature space gets truncated    
        

        case 'pca_denoise' %project onto top k PC to capture max variance, project back
            cutoffval = ops.preprocCutoff;
            [coeff,score,latent] = pca(fetj);
            cutoffvec = find(cumsum(latent)/sum(latent) > cutoffval);
            cutoff = cutoffvec(1);
            fet_proc(:,:,j) = score(:,1:cutoff)*coeff(:,1:cutoff)';

        otherwise

    end   
      
end


