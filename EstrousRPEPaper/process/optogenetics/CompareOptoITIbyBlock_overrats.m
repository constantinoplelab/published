function [nonopto_norm, opto_norm] =...
    CompareOptoITIbyBlock_overrats(ratlist, optoafterdates,...
    AlignedRats, SRats)

optosessions_rats = NaN(length(ratlist), 3);
nonoptosessions_rats = NaN(length(ratlist), 3);
for rat = 1:length(ratlist)

    ratname = ratlist{rat};
    disp([ratname ': ' num2str(rat) ' out of ' num2str(length(ratlist))])

    [optosessions_rats(rat, :), ~, nonoptosessions_rats(rat, :)] =...
        CompareOptoITIbyBlock(optoafterdates{rat},...
        AlignedRats.(ratname), SRats.(ratname));

end

%normalize to control sessions high block
nonopto_norm = nonoptosessions_rats./nonoptosessions_rats(:, 2);
opto_norm = optosessions_rats./nonoptosessions_rats(:, 2);

end