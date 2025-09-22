function [range_] = updateYlim(range, yl)
    
    if range(1)==0 && range(2)==0
        range_ = yl;
    else
        if yl(1) < range(1)
            range_(1) = yl(1);
        else
            range_(1) = range(1);
        end
        if yl(2) > range(2)
            range_(2) = yl(2);
        else
            range_(2) = range(2);
        end
    end
end