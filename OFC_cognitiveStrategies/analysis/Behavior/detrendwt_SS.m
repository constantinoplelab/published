function [A] = detrendwt_SS(A)

%initialize matrix for wait times
sess = nan(length(A.ntrials), max(A.ntrials)); 
ctr = 1;

wt = nan(length(A.wait_time), 1);
detrendedwt = wt;
wt(A.optout==1) = A.wait_time(A.optout==1);

xvec = (1:max(A.ntrials))';
s = []; x = [];
for jk = 1:length(A.ntrials)
    sess(jk,1:A.ntrials(jk)) = wt(ctr:ctr+A.ntrials(jk)-1);
    x = [x; xvec];
    s = [s; sess(jk,:)'];
    ctr = ctr+A.ntrials(jk);
end


%%we want to regress wait times against trial number, so we can subtract
%%it out later.

bad = find(isnan(s));
s(bad) = [];
x(bad) = [];
X = [x, ones(length(x),1)];
[b] = regress(s,X);
newy = (b(1).*xvec + b(2))'; 


ctr = 1;
for jk = 1:length(A.ntrials)
    detrendedwt(ctr:ctr+A.ntrials(jk)-1) = sess(jk,1:A.ntrials(jk))-newy(1:A.ntrials(jk));
    ctr = ctr+A.ntrials(jk);
end

A.wait_time = detrendedwt + b(2); 

end