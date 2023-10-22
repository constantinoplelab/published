function [pval, beta, goods, Y] = screen_perf(S)
%input: S struct, obtained from mysql.
%outputs:
 %p-value of slope coefficient for regression model relating wait time to 
 % reward volume, for each session.
 % beta is slope coefficient.
 % goods is a vector the sessions in which the animal's performance 
 % surpassed criterion. moving average
 % of slope must be positive, and the slope on single sessions must be 
 % positive for 2 consecutive days.
 % Y is moving average of slope coefficient in window of 7 sessions.
     %updated by cmc for readability, 10/19/21

pval = nan(length(S.pd),1);
beta = nan(length(pval),2);
wndw = 7; %moving window of 7 sessions for setting perf threshold

for j = 1:length(S.pd)
    if isstruct(S.pd{j})
        if sum(S.pd{j}.TrainingStage==9)==length(S.pd{j}.TrainingStage)...
                && length(S.pd{j}.TrainingStage)>100
            
            
            if isfield(S.pd{j}, 'wait_time')
                wt = S.pd{j}.wait_time;
                ctch = S.pd{j}.RewardDelay==100;
                rew = S.pd{j}.RewardAmount;
                wt(wt>mean(wt, 'omitnan')+2*std(wt, 'omitnan')) = nan;
                const = ones(length(wt),1);
                if sum(ctch==1)>5; %at least five catch trials
                    
                    %design matrix is reward volume and offset/constant term.
                    X = [rew(ctch==1), const(ctch==1)];
                    
                    %regress wait time against reward volume.
                    [beta(j,:),~,~,~,stats] = regress(wt(ctch==1),X);
                    pval(j) = stats(3);
                end     
                
            end
        end
    end
end

Y = movmean(beta(:,1), wndw, 'omitnan'); %moving average of slope

goods = find(Y>=0 & beta(:,1)>0); %moving average is >0 & that day has positive slope
d = [0;diff(goods)];
goods = goods(d==1); %must have at least 2 consecutive days of positive slope
