function plotPretty(x, y, color, errortype)
% Plot mean and shaded errorbar. Default error type is sem

sem = @(xx) nanstd(xx) ./ sqrt(size(xx,1));

if size(y,1) > 1
    if nargin<4
        shadedErrorBar(x, mean(y, 'omitnan'), sem(y), 'lineprops', ...
            {'color', color, 'linewidth', 1.25})
    elseif strcmpi(errortype, 'std') || strcmpi(errortype, 'stdev')
        shadedErrorBar(x, mean(y, 'omitnan'), nanstd(y), 'lineprops', ...
            {'color', color, 'linewidth', 1.25})
    end
elseif size(y,1) == 1 
    plot(x, y, 'color', color, 'LineWidth', 1.25); hold on
end

    set(gca, 'TickDir', 'out'); box off
    set(gca, 'FontSize', 13);
end