function [n_selective,p_selective] = selectiveNeuronCalc(A,mask1,mask2)
%SELECTIVENEURONCALC: decides via t test which neurons are selective
%for two different cases
%A: main data structure
%mask1/2: anonymous function handle that takes A and a neuron number to
%compute boolean of trials to consider for cases 1 and 2
%
% example: mask1 = @(A,j) A{j}.hits==1 gives all rewards on a trial
%
%output:
%n_selective: vector for selectivity. 0=not selective, 1=selective case 1,
%   -1 = selective case 2
%p_selective: vector of p values for the t-tests

nf = numel(A);
n_selective = zeros(nf,1);
p_selective = zeros(nf,1);

for j = 1:nf
    if A{j}.isUsable
        nspikes_1 = A{j}.nspikes(mask1(A,j));
        nspikes_2 = A{j}.nspikes(mask2(A,j));
        
        %test. using rate version
        %rates = mean(A{j}.hmat,2);
        %nspikes_1 = rates(mask1(A,j));
        %nspikes_2 = rates(mask2(A,j));
        
        [h,p] = ttest2(nspikes_1,nspikes_2);
        p_selective(j) = p;
        
        if h==1
            if nanmean(nspikes_1) > nanmean(nspikes_2)
                n_selective(j) = 1;
            else
                n_selective(j) = -1;
            end
        end       
    end
end