% run all of the preprocessing. save this for transparancy

epoch = 'reward';
decodertype = 'psth';
savename = strcat('~/projects/dynamics/data/maggie/preprocess_mlb_ofc_',epoch,'_',decodertype,'.mat')
run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)


epoch = 'reward';
decodertype = 'none';
savename = strcat('~/projects/dynamics/data/maggie/preprocess_mlb_ofc_',epoch,'_',decodertype,'.mat')
run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)


epoch = 'coff';
decodertype = 'psth';
savename = strcat('~/projects/dynamics/data/maggie/preprocess_mlb_ofc_',epoch,'_',decodertype,'.mat')
run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)


epoch = 'coff';
decodertype = 'none';
savename = strcat('~/projects/dynamics/data/maggie/preprocess_mlb_ofc_',epoch,'_',decodertype,'.mat')
run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)
%

epoch = 'son';
decodertype = 'psth';
savename = strcat('~/projects/dynamics/data/maggie/preprocess_mlb_ofc_',epoch,'_',decodertype,'.mat')
run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)


epoch = 'son';
decodertype = 'none';
savename = strcat('~/projects/dynamics/data/maggie/preprocess_mlb_ofc_',epoch,'_',decodertype,'.mat')
run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)

%% for psth 2

% run all of the preprocessing. save this for transparancy

epoch = 'reward';
decodertype = 'psth2';
savename = strcat('~/projects/dynamics/data/maggie/preprocess_mlb_ofc_',epoch,'_',decodertype,'.mat')
run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)


epoch = 'coff';
decodertype = 'psth2';
savename = strcat('~/projects/dynamics/data/maggie/preprocess_mlb_ofc_',epoch,'_',decodertype,'.mat')
run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)


epoch = 'son';
decodertype = 'psth2';
savename = strcat('~/projects/dynamics/data/maggie/preprocess_mlb_ofc_',epoch,'_',decodertype,'.mat')
run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)



