function standard_error_of_mean = sem(xx)

    standard_error_of_mean =...
        std(xx, 'omitnan')./sqrt(sum(~isnan(xx)));

end

