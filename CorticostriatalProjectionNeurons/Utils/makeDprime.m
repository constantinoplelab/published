function [d,dvar,dprime,h] = makeDprime(trials, goods, all_S,all_SU,all_index,event,sig,varargin)
% Inputs
%   trials: {:,2} cell array of trials for each good neuron
%   goods: subset of cells you want to use
%   all S
%   all SU
%   all_index
%   event
%   varargin: defaults to signed d', add value to make unsigned

%%output: 
%   d = difference between two psths
%   dvar = variance of each cell
%   dprime = dprime for each cell
%   h = 1 if there is a significant difference in mean frequency between the
%       two conditions. 

%keyboard
for ii = 1:length(goods)
session = all_index{goods(ii),4};
S = all_S{session};
SU = all_SU{goods(ii)};
T1 = trials{ii,1};
T2 = trials{ii,2};

T1Hmat = SU.hmat.(event)(T1,:);
T1PSTH(ii,:) = mean(T1Hmat,'omitnan');
n1 = length(T1);
T1Varaiance(ii,:) = var(T1Hmat,'omitnan');

T2Hmat = SU.hmat.(event)(T2,:);
T2PSTH(ii,:) = mean(T2Hmat,'omitnan');
n2 = length(T2);
T2Varaiance(ii,:) = var(T2Hmat,'omitnan');

%t-test on whether there is a significant difference in mean frequency
if sig==1
try
    h(ii,:) = ttest2(T1Hmat,T2Hmat);
catch
    h(ii,:) =NaN;
end
else
    h = NaN;
end

end %each cell

%difference in frequency on right-left
if ~isempty(varargin)
d = abs(T2PSTH-T1PSTH);
else
d = T2PSTH-T1PSTH;
end
dvar = sqrt((T1Varaiance+T2Varaiance)/2);
dprime = d./dvar;
