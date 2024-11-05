function [estradiol, RatConc, RatGroups, groups, m_mean, d_mean,...
    ep_mean, lp_mean, e_mean, m_err, d_err, ep_err, lp_err, e_err] =...
    serum_osmolality_estradiol(Outlier_factor, SerumOsmolalityTable,...
    EstradiolExp1Table, EstradiolExp2Table)

groups = {'male','diestrus','earlyproestrus','lateproestrus','estrus'};

%remove outliers
T = SerumOsmolalityTable;
T.Osmolality(T.Osmolality>mean(T.Osmolality)...
    +Outlier_factor*std(T.Osmolality)|T.Osmolality<...
    mean(T.Osmolality)-Outlier_factor*std(T.Osmolality)) = NaN; %removes two samples

%% As separate data points

%%%% Subset by group %%%%%
diestrus = T(strcmp(T.Group, 'diestrus'), :);
early_proestrus = T(strcmp(T.Group, 'earlyproestrus'), :);
late_proestrus = T(strcmp(T.Group, 'lateproestrus'), :); %3PM or later
estrus = T(strcmp(T.Group, 'estrus'), :);
male = T(strcmp(T.Group, 'male'), :);

%%Get averages and SEM
d_mean = mean(diestrus.Osmolality, 'omitnan');
ep_mean = mean(early_proestrus.Osmolality, 'omitnan');
lp_mean = mean(late_proestrus.Osmolality, 'omitnan');
e_mean = mean(estrus.Osmolality, 'omitnan');
m_mean = mean(male.Osmolality, 'omitnan');

d_err = std(diestrus.Osmolality, 'omitnan')./sqrt(sum(~isnan(diestrus.Osmolality)));
ep_err = std(early_proestrus.Osmolality, 'omitnan')./sqrt(sum(~isnan(early_proestrus.Osmolality)));
lp_err = std(late_proestrus.Osmolality, 'omitnan')./sqrt(sum(~isnan(late_proestrus.Osmolality)));
e_err = std(estrus.Osmolality, 'omitnan')./sqrt(sum(~isnan(estrus.Osmolality)));
m_err = std(male.Osmolality, 'omitnan')./sqrt(sum(~isnan(male.Osmolality)));

%% Avg by rat
%%Get averages per rat
ratlist = unique(T.RatID);
RatNames = cell(length(ratlist), 1);
RatGroups = cell(length(ratlist), 1);
RatConc = NaN(length(ratlist), 1);
for rat = 1:length(ratlist)
    ratT = T(strcmp(T.RatID, ratlist(rat)), :);
    RatNames(rat) = ratlist(rat);
    RatGroups(rat) = ratT.Group(1);
    RatConc(rat) = mean(ratT.Osmolality, 'omitnan');
end

estradiol = NaN(length(ratlist), 1);
estr_exp1 = EstradiolExp1Table;
estr_exp2 = EstradiolExp2Table;
for rat = 1:length(ratlist)
    ratname = ratlist(rat);
    if ismember(ratname, estr_exp1.Rat) && ismember(ratname, estr_exp2.Rat)
        e1 = estr_exp1.Concentration(ismember(estr_exp1.Rat, ratname));
        e2 = estr_exp2.Concentration(ismember(estr_exp2.Rat, ratname));
        estradiol(rat) = mean([e1 e2]);
    elseif ismember(ratname, estr_exp1.Rat)
        estradiol(rat) = estr_exp1.Concentration(ismember(estr_exp1.Rat, ratname));
    elseif ismember(ratname, estr_exp2.Rat)
        estradiol(rat) = estr_exp2.Concentration(ismember(estr_exp2.Rat, ratname));
    end
end

%292395's E2 expression is an outlier as a rat in estrus
%remove it
estradiol(7) = NaN; %<1 std below the mean of all samples

end