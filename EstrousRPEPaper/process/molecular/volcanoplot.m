function [log2FC, y, txt1, txt2, txt3, stagescompared, stageslfc] =...
    volcanoplot(MS, comparison, direction)
%input: comparison = 1 for diestrus vs proestrus
                    %3 for diestrus vs estrus

if comparison == 1
    stagescompared = 'Diestrus vs. Proestrus';
    stageslfc = 'Diestrus/Proestrus';
    txt1 = 'Slc6a3';
    txt2 = 'Slc6a4';
    txt3 = 'Smpd3';
    log2FCcolumn = 12;
    if direction == -1
        stagescompared = 'Proestrus vs. Diestrus';
        stageslfc = 'Proestrus/Diestrus';
    end
elseif comparison == 2
    stagescompared = 'Proestrus vs. Estrus';
    stageslfc = 'Proestrus/Estrus';
    txt1 = 'Slc6a4';
    txt2 = 'Smpd3';   
    if direction == -1
        stagescompared = 'Estrus vs. Proestrus';
        stageslfc = 'Estrus/Proestrus';
    end
elseif comparison == 3
    stagescompared = 'Diestrus vs. Estrus';
    stageslfc = 'Diestrus/Estrus';
    txt1 = 'Slc6a4';
    txt2 = 'Smpd3';
    txt3 = 'Slc6a3';  
    log2FCcolumn = 11;
    if direction == -1
        stagescompared = 'Estrus vs. Diestrus';
        stageslfc = 'Estrus/Diestrus';
    end    
end
                    
%define Y (p-values)
y = MS.x_Log10P_value;

%define X (log fold change)
log2FC = MS{:, log2FCcolumn};
if direction == -1
    log2FC = log2FC*-1; %if you want to change the direction of the comparison
end

end