function meas = fun_clusterConsistency(lab1,lab2)
%function to check cluster consistency, relative to labeling 2, as
%MI(lab1,lab2)/H(lab2)

%sample script to show how the similarity measure changes with probability
unique1 = unique(lab1);
unique2 = unique(lab2);
n1 = numel(unique1); %classes in lab1
n2 = numel(unique2); %classes in lab2
ns = numel(lab1); %number samples

p2 = zeros(n2,1);
p1 = zeros(n1,1);
pcond = zeros(n2,n1); %conditional probabilities, p(L2|L1)

%marginal probabilities
for j = 1:n2
    p2(j) = sum(lab2==unique2(j))/ns;
end
for j = 1:n1
    p1(j) = sum(lab1==unique1(j))/ns;
end

%conditional probabilities
for j = 1:n2
    for k = 1:n1
       pcond(j,k) = sum(lab2==unique2(j) & lab1==unique1(k))/sum(lab1==unique1(k));        
    end
end

%deal with missing cases and inf/nan. make very small
pcond(pcond==0) = 1e-6;

HC =  -p2'*log2(p2);
Hcond = -sum(pcond.*log2(pcond))*p1;
%Hcond = -pC(1)*pcond(1,:)*log2(pcond(1,:)') -pC(2)*pcond(2,:)*log2(pcond(2,:)');

MI = HC-Hcond;

meas = MI/HC;




 