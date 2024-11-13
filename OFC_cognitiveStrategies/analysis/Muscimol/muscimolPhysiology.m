function out = muscimolPhysiology(muscimolPhysDataPath)

sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));

SU_L = load([muscimolPhysDataPath, filesep, 'M011_L.mat']);
SU_R = load([muscimolPhysDataPath, filesep, 'M011_R.mat']);
T = readtable([muscimolPhysDataPath, filesep, 'acuteMuscimol_info.xlsx']);
left = strcmpi(T.Rat, 'M011_L');
right = strcmpi(T.Rat, 'M011_R');

pre_time = T.pre_time(1) *60; %baseline recording before local muscimol infusion
post_time = T.post_time(1) *60; %recording time after muscimol infusion
dividePoint = pre_time : pre_time + (T.infusion_time(1)*60); %total recording time

% split spikes into pre and post muscimol -- Left recording
mat_L = nan(length(SU_L.SU),3);

for u = 1:length(SU_L.SU)
    %calculate distance from center of muscimol infusion to channel
    mat_L(u,1) = sqrt((T.distance(left) * 1000)^2 + SU_L.SU{u}.channel_depth^2); 
    
    %caluculate pre- and post-muscimol firing rates
    mat_L(u,2) = numel(find(SU_L.SU{u}.st < dividePoint(1))) / pre_time; % pre-muscimol firing rate 
    mat_L(u,3) = numel(find(SU_L.SU{u}.st > dividePoint(end))) / post_time; % post-muscimol firing rate
end

% split spikes into pre and post muscimol -- Right recording
mat_R = nan(length(SU_R),3);

for u = 1:length(SU_R.SU)
    %calculate distance from center of muscimol infusion to channel
    mat_R(u,1) = sqrt((T.distance(right) * 1000)^2 + SU_R.SU{u}.channel_depth^2); 
    
    %caluculate pre- and post-muscimol firing rates
    mat_R(u,2) = numel(find(SU_R.SU{u}.st < dividePoint(1))) / pre_time; % pre-muscimol firing rate 
    mat_R(u,3) = numel(find(SU_R.SU{u}.st > dividePoint(end))) / ...
        (post_time-T.infusion_time(1)); % post-muscimol firing rate
end

mat = [mat_L; mat_R];

mat(mat(:,2) < 0.4, :) = []; %Exclude cells with a baseline firing rate <0.4 hz

[~,idx] = sort(mat(:,1)); 
sortedmat = mat(idx,:);
sortedmat(sortedmat(:,2) == 0) = eps;
sortedmat(sortedmat(:,3) == 0) = eps;

ratio = sortedmat(:,3) ./ sortedmat(:,2); 

bins = sortedmat(1,1):100:sortedmat(end,1); %original was 100
out.binnedratio = arrayfun(@(x) mean(ratio(find(sortedmat(:,1) >= bins(x) & ...
    sortedmat(:,1) < bins(x+1))), 'omitnan'), 1:length(bins)-1);
out.binnedEr = arrayfun(@(x) sem(ratio(find(sortedmat(:,1) >= bins(x) & ...
    sortedmat(:,1) < bins(x+1)))), 1:length(bins) -1);

%fit sigmoid
remove = isnan(out.binnedratio);
out.bins = bins(1:end-1)./1000;
out.x = out.bins; out.x(remove) = [];
y = out.binnedratio; y(remove) = [];
[finalParams, MSE, initCond, fitParams] = my_fit_sigmoid(out.x, y, 10);

out.sigFit = finalParams(4) + (finalParams(1)-finalParams(4)) ./ ...
    (1 + exp(-finalParams(2)*(out.x-finalParams(3))));

% example cell
egCell = SU_R.SU{5};
totalTime = T.pre_time(1)*60 + T.post_time(1)*60;

wndw = [0 totalTime];
bins = 60; %1 minute bins

xvec = [wndw(1):bins:wndw(2)];
d = diff(xvec)/2;
edges = [xvec(1)-d(1), xvec(1:end-1)+d, xvec(end)+d(end)];

[n, ~] = histcounts(egCell.st, edges);   
hmat = n./bins;
out.hmat = smooth(hmat);
out.xvec = (-pre_time:bins:post_time)'./60;


