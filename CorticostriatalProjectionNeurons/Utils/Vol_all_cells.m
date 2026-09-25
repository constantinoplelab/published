function [avgSpk] = Vol_all_cells(all_SU,goods,all_S,all_index,win,varargin)
%plot heatmats aligned to reward volume for all cells
%plots the average psth for reward volume at each task event. 
%plots a hmat

% input:
% all_SU - the data of all cells
% goods - cells you want to include
% all_S - data of all cells
% all_index - indexing of all cells
% win - [start,end] of window you want plotted
% varargin:
    % alignto to sort the data too, otherwise defaults to time to peak ff on COff
    % 1 = z-score trial by trial
    % define the index by which to sort a hmat of volume encoding
    % plot FR by sorted index
% output:
% figure of heatmat with z-scored ff in aligned to reward volume

load("\\constantinoplelab.cns.nyu.edu\server\Maggie\Software\Matlab\Colormaps\TriColor.mat")
vcolors = [.6 .2 .9; .1 .5 .9; 0 .4 0; 1 .65 0; 1 0 0]; %need to play with colors to be color blind friendly

%% z-score each cell, each trial
nanzscore = @(x)(x - mean(x,2, 'omitnan'))./std(x,0,2, 'omitnan');
alignto = {'COFF', 'SON','SOFF','Rew','Opt','CON'};



Label_vol = {'5uL','10uL','20uL','40uL','80uL'};

for a = 1:length(goods)
        clusternames(a) = all_SU{1,goods(a)}.cluster_id;
end

for k = 1:5
    if nargin>6  && varargin{2}==1
    % for i = 1:length(goods)
    %     all_SU{1,goods(i)}.hmat.(string(alignto(k))) = nanzscore(all_SU{1,goods(i)}.hmat.(string(alignto(k))));
    % end
    end
    %% calculate trial-average of z-scored spike rate for each cell
    for i = 1:length(goods)
        this =  all_index{goods(i),4};  %S struct number
        %trials = all_behEvents{1,this}.NoVio;
        %% find trials for each reward vol
        S = all_S{1,this};
        [S.RewardAmount] = convertreward(S.RewardAmount);

        %% find volume response
        for v = 1:5
            %vol, non-violation, mixed block only
             V = S.vios==0 & S.Block==1 & S.RewardAmount==v;

            % %vol, non-violation, all blocks
            % V = S.vios==0 & S.RewardAmount==v;
            
             % %vol, violation, all blocks
             %V = S.RewardAmount==v;

             % %vol, all blocks, only hits
             %V = S.RewardAmount==v & S.hits==1;

             % %vol, optouts, all blocks
             %V = S.RewardAmount==v & S.vios==0 & S.optout==1;

            Vol{v,:} = find(V==1);

            avgSpk.(string(alignto{k})){v}(i,:) = mean(all_SU{1,goods(i)}.hmat.(string(alignto(k)))(Vol{v},:),1, 'omitnan'); % average spike rate for a cell over all trials in a session

        end
    end
end

if nargin>9
else
%% plotting
%T = all_SU{1,1}.xvec;
figure('PaperSize',[11 8.5],'color','w')
xvec = all_SU{goods(1)}.xvec.CON;
T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');


%sorting
if isempty(varargin)
    vol5 = avgSpk.(alignto{1}){1, 1};
else
    vol5 = avgSpk.(varargin{1}){1, 1};
end

vol5_top = max(vol5,[],2);
for i = 1:length(vol5_top)
    vol5_local(i)= find(vol5(i,:)==vol5_top(i),1);
end
vol5_local_sort = [vol5_local;1:length(vol5_top)]';
vol5_local_sort = sortrows(vol5_local_sort);
vol5_indx = vol5_local_sort(:,2);

for k=1:5
    for v = 1:5
        p(k) = subplot(5,5,(5*(v-1)+k));
        imagesc(xvec(T), linspace(1,size(goods,2)), avgSpk.(string(alignto{k})){v}(vol5_indx,T))

        if v==1
            title(alignto(k))
        end

        if nargin>6 && varargin{2}==1
        clim([-2 2])
        else
        clim([0 30])
        end
        xline(0,'--w',LineWidth=1)
        box off
        if k==5
            colorbar
            p(k).Position = [(p(k).Position(1)) (p(k).Position(2)) (p(k-1).Position(3)) (p(k-1).Position(4))];
            drawnow
            yticklabels([])
        elseif k==1
            yticklabels(clusternames)
            ylabel(Label_vol{v})
        else
            yticks([])
            yticklabels([])
        end
    end
    if nargin>6 && varargin{2}==1
    colormap(TriColor)
    end
end
    


if nargin>7
    if varargin{3}==1
    %zscore after averaging
    for k = 1:5
        for i = 1:length(goods)
            for v = 1:5
                hmat(v,:) = avgSpk.(alignto{k}){v}(i,:);
                hmatm = mean(hmat,'omitnan');
                hmatstd = std(hmat,'omitnan');
                nanzscore2 = @(x)(x - hmatm./hmatstd);
            end
            zspk = nanzscore2(hmat);

            % figure('color','w')
            % for v = 1:5
            %     subplot(1,2,1)
            % plot(zspk(v,T),'Color',vcolors(v,:))
            % hold on
            %     subplot(1,2,2)
            % plot(hmat(v,T),'Color',vcolors(v,:))
            % hold on
            % end

            for v = 1:5
                avgSpk.(alignto{k}){v}(i,:) = zspk(v,:);
            end
        end
    end
    end
end


limit = 0;

%% Plot mean response
figure('color','w')
for k = 1:5
    subplot(1,5,k)

    for vol = 1:5
        % avgSpk.(alignto{k}){vol} = nanzscore(avgSpk.(alignto{k}){vol});
        plot(xvec(T),mean(avgSpk.(alignto{k}){vol}(:,T),'omitnan'),'Color',vcolors(vol,:))
        limit = [limit,mean(avgSpk.(alignto{k}){vol}(:,T),'omitnan')];
        % shadedErrorBar(xvec(T),mean(avgSpk.(alignto{k}){vol}(:,T),'omitnan'),sem(avgSpk.(alignto{k}){vol}(:,T),'omitnan'),'lineprops',{'Color',vcolors(vol,:)})
        hold on;
    end
end
subplot(1,5,1);
if nargin>6  && varargin{2}==1
    ylabel('Firing rate (z-scored)')
else
    ylabel('Firing rate (Hz)')
end
for k = 1:5
    subplot(1,5,k)
    top = max(limit(1,2:end));
    bottom = min(limit(1,2:end));
    ylim([bottom-abs(bottom/10) top+top/50])
    title(alignto{k})
    box off; xlabel('Time (s)');set(gca,'TickDir','out');xline(0,'--');
end
legend(Label_vol,'location','northeast')


%% Plot FR for mixblock
if nargin>8
ind = varargin{4};
t = {'5uL','10uL','20uL','40uL','80uL'};
figure('color','w')
tiledlayout(1,5)
%all tiles
for a = 1:5
    nexttile
imagesc(xvec(T),1:length(goods),avgSpk.Rew{a}(ind,T))
colorbar;
if nargin>6 && varargin{2}==1
clim([-1 3])
else
clim([0 30]); 
end
xline(0,'--w');xlabel('Time to Reward (s)')

title(t{a});box off;set(gca,'TickDir','out')
end
%tile one
nexttile(1)
yticks(1:length(goods));yticklabels(clusternames(ind));ylabel('Cells');

%rest of tiles
for k = 2:5
    nexttile(k)
    yticks(1:length(goods));yticklabels([]);
end

sgtitle('Mixed block volume encoding at Reward')
end

% %% 
% keyboard
% figure
% for i = 1:length(goods)
%     subplot(4,2,i)
%     for a = 1:5
%         plot(xvec(T),avgSpk.Rew{a}(i,T),'color',vcolors(a,:));hold on;
%     end
%     box off;
% end
end
