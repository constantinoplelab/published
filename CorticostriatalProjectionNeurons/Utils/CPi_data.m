function Avg = CPi_data
%generates the CPi figure data


%% load Jonny DCS
load("Z:\histology\Jonny\Analysis\RUN25 heatmats\RUN25_slide1 heatmats\RUN25_1_3.mat")
R1 = Count.VTA;

load("Z:\histology\Jonny\Analysis\RUN26 heatmats\RUN26_slide1 heatmats\RUN26_1_2.mat")
R2 = Count.VTA;

load("Z:\histology\Jonny\Analysis\X012 heatmats\X012_slide1 heatmats\X012_1_1.mat")
R3 = Count.VTA;

%% load Maggie/Mitzi DLS
load('Z:\histology\Maggie\Retrograde_striatum_CTB\Protocol Information\Combined_data\2022-11-07.mat');

location ={'DLS'};
i = 1;
    
    % G035
    a = cat(4,C.G035.slide_2.section_B.(string(location(i))), C.G035.slide_2.section_D.(string(location(i))), ...
        C.G035.slide_2.section_E.(string(location(i))));
    R4 = mean(a,4,'omitnan');
    
    % S020
    a = cat(4, C.S020.slide_2.section_A5.(string(location(i))), C.S020.slide_2.section_A5p.(string(location(i))));
    R5 = mean(a,4,'omitnan');
    
    
    % RUN09
    a = cat(4, C.RUN09.slide_3.section_A_01.(string(site(i))), C.RUN09.slide_3.section_A_02.(string(site(i))));
    C.RUN09.slide_3.average.A.(string(site(i))) = mean(a,4,'omitnan');
    a = cat(4, C.RUN09.slide_3.section_G_01.(string(site(i))), C.RUN09.slide_3.section_G_02.(string(site(i))));
    C.RUN09.slide_3.average.G.(string(site(i))) = mean(a,4,'omitnan');
    a = cat(4, C.RUN09.slide_3.section_I_01.(string(site(i))), C.RUN09.slide_3.section_I_02.(string(site(i))));
    C.RUN09.slide_3.average.I.(string(site(i))) = mean(a,4,'omitnan');
    a = cat(4, C.RUN09.slide_3.section_J_01A.(string(site(i))), C.RUN09.slide_3.section_J_01B.(string(site(i))), ...
        C.RUN09.slide_3.section_J_02_A.(string(site(i))), C.RUN09.slide_3.section_J_02_B.(string(site(i))));
    C.RUN09.slide_3.average.J.(string(site(i))) = mean(a,4,'omitnan');
    a = cat(4, C.RUN09.slide_3.section_K_01.(string(site(i))), C.RUN09.slide_3.section_K_02.(string(site(i))));
    C.RUN09.slide_3.average.K.(string(site(i))) = mean(a,4,'omitnan');
    
    a = cat(4, C.RUN09.slide_3.average.A.(string(site(i))), C.RUN09.slide_3.section_D.(string(site(i))), ...
        C.RUN09.slide_3.section_E.(string(site(i))), C.RUN09.slide_3.section_H.(string(site(i))));
    R6 = mean(a,4,'omitnan');
    
    


R = cat(3,R1,R2,R3,R4,R5,R6);
Avg = mean(R,3,'omitnan');
Ro = cat(3,R1,R2,R3);
Avgo = mean(Ro,3,'omitnan');

load('Z:\histology\Maggie\Retrograde_striatum_CTB\Protocol Information\Combined_data\colormap.mat');
% figure
% ax{1}=subplot(1,2,1);
% imagesc(Avg)
% ax{2}=subplot(1,2,2);
% imagesc(Avgo)
% for a = 1:2
% set(ax{a}, 'Colormap', cmap);
% end




