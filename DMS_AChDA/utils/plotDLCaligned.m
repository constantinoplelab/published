function plotDLCaligned(AlignmentStruct)

figure;
set(gcf, 'units', 'inches', 'position', [4.5,5.77,10.66,3.406])
A = fields(AlignmentStruct.speed);
range1 = [0 0];
for a=1:length(A)
    % head speed
    subplot(2, length(A), a)
    plotPretty(AlignmentStruct.Times, AlignmentStruct.speed.(A{a}), ...
        'k', 'std')
    if a==1
        ylabel('Head speed (cm/s)')
    end
    yl = ylim;
    range1 = updateYlim(range1, yl);
    title(A{a})
    
    subplot(2, length(A), length(A)+a)
    plotPretty(AlignmentStruct.Times, AlignmentStruct.headingCenter.(A{a}), ...
        'k', 'std')
    set(gca, tickdir='in')
    xlim([-0.5 1])
    ylim([0 180])
    if a==1
        ylabel('Heading w.r.t. CP')
    end
    axis square
    xline(0, 'k--')
end

for a=1:length(A)
    subplot(2, length(A), a)
    ylim(range1)
    xlim([-0.5 1])
    axis square
    xline(0, 'k--')
end

end