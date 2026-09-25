
% the commands to make the preprocess data

usetorch = false;
savedir = '~/projects/dynamics/data/maggie/';

%% reward epoch
disp('reward epochs')
epoch = 'reward';
decodertype = 'psth';
usemlb = true;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

epoch = 'reward';
decodertype = 'none';
usemlb = true;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)


epoch = 'reward';
decodertype = 'psth';
usemlb = false;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

epoch = 'reward';
decodertype = 'none';
usemlb = false;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

% coff
disp('coff epochs')
epoch = 'coff';
decodertype = 'psth';
usemlb = true;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

epoch = 'coff';
decodertype = 'none';
usemlb = true;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)


epoch = 'coff';
decodertype = 'psth';
usemlb = false;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

epoch = 'coff';
decodertype = 'none';
usemlb = false;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

% son
disp('Son epochs')
epoch = 'son';
decodertype = 'psth';
usemlb = true;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

epoch = 'son';
decodertype = 'none';
usemlb = true;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)


epoch = 'son';
decodertype = 'psth';
usemlb = false;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

epoch = 'son';
decodertype = 'none';
usemlb = false;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

%% the psth 2 ones

disp('reward epochs')
epoch = 'reward';
decodertype = 'psth2';
usemlb = true;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)


epoch = 'reward';
decodertype = 'psth2';
usemlb = false;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

% coff
disp('coff epochs')
epoch = 'coff';
decodertype = 'psth2';
usemlb = true;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

epoch = 'coff';
decodertype = 'psth2';
usemlb = false;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

% son
disp('Son epochs')
epoch = 'son';
decodertype = 'psth2';
usemlb = true;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

epoch = 'son';
decodertype = 'psth2';
usemlb = false;
if usemlb
    savename = strcat(savedir,'preprocess_mlb_shuffle_',epoch,'_',decodertype,'.mat');
else
    savename = strcat(savedir,'preprocess_mlb2_shuffle_',epoch,'_',decodertype,'.mat');
end
disp(savename)
run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

