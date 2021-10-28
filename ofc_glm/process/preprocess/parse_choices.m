function [chosenprob, chosenval, hits, notchosenprob, notchosenval] = parse_choices(S)
%PARSE_CHOICES determines if trial was rewarded, and probability  and volume  of reward
%
%inputs: S data structure with num flashes, num beeps, hits, etc.
%
%outputs:   chosenprob: probability of reward in chosen port
%           chosenval: volume of water reward in chosen port (mL)
%           hits: 0 = not rewarded, 1 = rewarded, nan = violated 
%                (aka did not particpate in trial and left center poke)
%           notchosenprob: prob. of reward from non-chosen port
%           notchosenval: volume of water reward from non-chosen port (mL)
%
%output shapes will all be (n x 1) arrays, where n is number of trials in 
%a session

%initialize
chosenprob = nan(length(S.pd{1}.hits), 1);
chosenval = chosenprob;
notchosenprob = chosenprob;
notchosenval = chosenprob;

%find values or probability and water volume, case of rat choosing right port
r = find(S.pd{1}.went_right==1);
chosenprob(r) = S.pd{1}.right_prob(r);
chosenval(r) = S.pd{1}.this_right_volume(r);
notchosenprob(r) = S.pd{1}.left_prob(r);
notchosenval(r) = S.pd{1}.this_left_volume(r);

%similarly for rat choosing left port
l = find(S.pd{1}.went_right==0);
chosenprob(l) = S.pd{1}.left_prob(l);
chosenval(l) = S.pd{1}.this_left_volume(l);
notchosenprob(l) = S.pd{1}.right_prob(l);
notchosenval(l) = S.pd{1}.this_right_volume(l);

%determine if trial was rewarded by matching port choice with a port hit
hits = nan(length(chosenprob), 1);
h = find(S.pd{1}.went_right==1 & S.pd{1}.right_hit==1 | ...
    S.pd{1}.went_left==1 & S.pd{1}.left_hit==1);
hits(h) = 1;
h = find(S.pd{1}.went_right==1 & S.pd{1}.right_hit==0 | ...
    S.pd{1}.went_left==1 & S.pd{1}.left_hit==0);
hits(h) = 0;