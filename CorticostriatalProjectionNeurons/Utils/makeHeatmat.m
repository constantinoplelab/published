function SU = makeHeatmat(SU, S, behEvents, wndw, bins, varargin)
% input
% SU struct
% S struct
% behEvents struct
% wndw = [-2,4] usually
% bin = 0.05s usually
% varargin = include raster or not

if nargin<5
    bins = 0.05;
end
%side LED on
sideEvents = {S.Lled(:,1) S.Rled(:,1) S.l_opt(:,1) S.r_opt(:,1)};   %pull times of all events
[s,~] = find(cell2mat(arrayfun(@(x) ~isnan(sideEvents{x}), 1:4, 'uniformoutput', false)));  %remove all nans
s = sort(s);
s_times = [S.Lled(:,1); S.Rled(:,1);S.l_opt(:,1);...
    S.r_opt(:,1)];
s_times= sort(s_times(~isnan(s_times)));
% s_ON = nan(length(S.Cled),1);
% s_ON(s) = s_times;

%reward
rewEvents = {S.Lled(:,3),S.Rled(:,3)};
[r, ~] = find(cell2mat(arrayfun(@(x) ~isnan(rewEvents{x}), 1:2, 'uniformoutput', false)));
r = sort(r);
r_times =  [S.Lled(:,3),S.Rled(:,3)];
r_times = sort(r_times(~isnan(r_times)));
rew = nan(length(S.Cled),1);
rew(r) = r_times;

%side LED off
soff_times = [S.Lled(:,2); S.Rled(:,2);]; %only use side off for rewarded trials
soff_times = sort(soff_times(~isnan(soff_times)));
s_OFF = nan(length(S.Cled),1);
s_OFF(r) = soff_times;

%opt-out
optEvents = {S.l_opt(:,3),S.r_opt(:,3)};
[o, ~] = find(cell2mat(arrayfun(@(x) ~isnan(optEvents{x}), 1:2, 'uniformoutput', false)));
o = sort(o);
o_times = [S.l_opt(:,3),S.r_opt(:,3)];
o_times = sort(o_times(~isnan(o_times)));
opt = nan(length(S.Cled),1);
opt(o) = o_times;

%CLED on, CLED off, side LED on, side LED off, reward, opt out
align0 = {S.Cled(:,1),S.Cled(:,2),s_times,soff_times,r_times,o_times};
cond = {behEvents.All,behEvents.All, s, behEvents.Rewarded, ...
    behEvents.Rewarded, behEvents.Optout};

for ii = 1:length(SU)
    for c = 1:length(align0)
        if isempty(varargin)

            [xvec,hmat, ~] = alignSpikes(SU{ii}.st, cond{c}, wndw, bins, ...
                align0{c}, S.Cled);
            if c == 1
                SU{ii}.hmat.CON = hmat; %CLED on hmat
                SU{ii}.xvec.CON = xvec;
            elseif c == 2
                SU{ii}.hmat.COFF = hmat; %CLED off hmat
                SU{ii}.xvec.COFF = xvec;
                SU{ii}.raster.COFF = cell(length(S.Cled),1);
            elseif c == 3
                SU{ii}.hmat.SON = nan(length(S.Cled),length(xvec));
                SU{ii}.hmat.SON(s,:) = hmat; %side LED on
                SU{ii}.xvec.SON = xvec;
            elseif c == 4
                SU{ii}.hmat.SOFF = nan(length(S.Cled), length(xvec));
                SU{ii}.hmat.SOFF(r,:) = hmat;
                SU{ii}.xvec.SOFF = xvec;
            elseif c == 5
                SU{ii}.hmat.Rew = nan(length(S.Cled),length(xvec));
                SU{ii}.hmat.Rew(r,:) = hmat; %Reward
                SU{ii}.xvec.Rew = xvec;
            elseif c == 6
                SU{ii}.hmat.Opt = nan(length(S.Cled),length(xvec));
                SU{ii}.hmat.Opt(o,:) = hmat; %Opt out
                SU{ii}.xvec.Opt = xvec;
            end
        else
            
            [xvec,hmat, raster] = alignSpikes(SU{ii}.st, cond{c}, wndw, bins, ...
                align0{c}, S.Cled);
            if c == 1
                SU{ii}.hmat.CON = hmat; %CLED on hmat
                SU{ii}.xvec.CON = xvec;
                SU{ii}.raster.CON = raster;
            elseif c == 2
                SU{ii}.hmat.COFF = hmat; %CLED off hmat
                SU{ii}.xvec.COFF = xvec;
                SU{ii}.raster.COFF = cell(length(S.Cled),1);
                SU{ii}.raster.COFF = raster;
            elseif c == 3
                SU{ii}.hmat.SON = nan(length(S.Cled),length(xvec));
                SU{ii}.hmat.SON(s,:) = hmat; %side LED on
                SU{ii}.xvec.SON = xvec;
                SU{ii}.raster.SON = cell(length(S.Cled),1);
                SU{ii}.raster.SON(s,:) = raster;
            elseif c == 4
                SU{ii}.hmat.SOFF = nan(length(S.Cled), length(xvec));
                SU{ii}.hmat.SOFF(r,:) = hmat;
                SU{ii}.xvec.SOFF = xvec;
                SU{ii}.raster.SOFF = cell(length(S.Cled),1);
                SU{ii}.raster.SOFF(r,:) = raster;
            elseif c == 5
                SU{ii}.hmat.Rew = nan(length(S.Cled),length(xvec));
                SU{ii}.hmat.Rew(r,:) = hmat; %Reward
                SU{ii}.xvec.Rew = xvec;
                SU{ii}.raster.Rew = cell(length(S.Cled),1);
                SU{ii}.raster.Rew(r,:) = raster;
            elseif c == 6
                SU{ii}.hmat.Opt = nan(length(S.Cled),length(xvec));
                SU{ii}.hmat.Opt(o,:) = hmat; %Opt out
                SU{ii}.xvec.Opt = xvec;
                SU{ii}.raster.Opt = cell(length(S.Cled),1);
                SU{ii}.raster.Opt(o,:) = raster;
            end
        end

    end
end
