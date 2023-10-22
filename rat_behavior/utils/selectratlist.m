a = load('Z:\ProcessedRatData\ratList.mat');
ratListAll = a.ratList;
sexListAll = a.sexList;

ratList = {};
sexList = {};

test = nan(size(ratListAll));

savedir = 'Z:\ProcessedRatData\A_Structs_Final\';
for rr = 1:length(ratListAll)
    disp(rr)
    a = load(['Z:\ProcessedRatData\A_Structs\ratTrial_' ratListAll{rr}]);
    A = a.A;

    test(rr) = all(A.test_block==40 & A.adapt_block==40);

    if all(A.test_block==40 & A.adapt_block==40)
        ratList{end+1} = ratListAll{rr}; %#ok<*SAGROW>
        sexList{end+1} = sexListAll{rr};
        save([savedir 'ratTrial_' ratListAll{rr} '.mat'], 'A')
     end

end
ratList = ratList';

save([savedir 'ratList.mat'], 'ratList', 'sexList') 
