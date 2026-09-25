function p = compareDeltaWT_CPiCPr

datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure1\Data';

%load data OFC-CPi
CPi = load(fullfile(datpath,'CPi_opto.mat'));
[highc_CPi,highe_CPi,lowc_CPi,lowe_CPi,mixc_CPi,mixe_CPi]= wt_notsided(CPi.ExpA,CPi.ConA,CPi.Raw_data);
%change
lowchange_CPi = lowe_CPi-lowc_CPi;
highchange_CPi = highe_CPi-highc_CPi;
mixchange_CPi = mixe_CPi-mixc_CPi;
%reshape
lowchange_CPi = reshape(lowchange_CPi,[55,1]);
highchange_CPi = reshape(highchange_CPi,[55,1]);
mixchange_CPi = reshape(mixchange_CPi,[55,1]);
deltaCPi = [lowchange_CPi;mixchange_CPi;highchange_CPi];

%load data OFC-CPr
% [ExpA_CPr,ConA_CPr,Raw_data_CPr] = Opto_dailycheck('VS',1,0);
CPr = load(fullfile(datpath,'CPr_opto.mat'));
[highc_CPr,highe_CPr,lowc_CPr,lowe_CPr,mixc_CPr,mixe_CPr]= wt_notsided(CPr.ExpA,CPr.ConA,CPr.Raw_data);
%change
lowchange_CPr = lowe_CPr-lowc_CPr;
highchange_CPr = highe_CPr-highc_CPr;
mixchange_CPr = mixe_CPr-mixc_CPr;
%reshape
lowchange_CPr = reshape(lowchange_CPr,[35,1]);
highchange_CPr = reshape(highchange_CPr,[35,1]);
mixchange_CPr = reshape(mixchange_CPr,[35,1]);
deltaCPr = [lowchange_CPr;mixchange_CPr;highchange_CPr];

%Compare the delta WT for CPi and CPr neurons
p = ranksum(deltaCPi,deltaCPr);
