function ind = plotdprime(hmat,T,xvec,varargin)
% inputs:
    % hmat: hmat of the data you want to plot
    % T: index of hmat you want to plot
    % times: axes times corresponding to the T values of the hmat
    % varagrin =
    % 1 sort time to peak
    % 1 win1 2 win2

load("\\constantinoplelab.cns.nyu.edu\server\Maggie\Software\Matlab\Colormaps\TriColor.mat")
%times = round(times,2);

figure('color','w')
if nargin<4
imagesc(xvec(T),1:length(hmat(:,1)),hmat(:,T))
title('dprime by cell')

elseif nargin==4
[a,i]= min(hmat(:,T),[],2);
[a,ind] = sort(i,'descend');
imagesc(xvec(T),1:length(hmat(:,1)),hmat(ind,T))

elseif nargin>4
% sort it
v0 = hmat(:,xvec>=varargin{1} & xvec<=varargin{2});
v1 = find(v0==-Inf | v0==Inf);
v0(v1) = NaN;
v0 = mean(v0',1,'omitnan')';
[ind2,ind]=sort(v0,'descend');
imagesc(xvec(T),1:length(hmat(:,1)),hmat(ind,T))
title('dprime by cell sorted by peak')

end
yticks([]);xlabel('time (s)');ylabel('cells');xline(0,'--');box off; set(gca,'tickdir','out');
colorbar
clim([-1 1])
colormap(TriColor)


