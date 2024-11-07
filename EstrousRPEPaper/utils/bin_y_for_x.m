function [ymean_binned, x_binned, yer_binned] = bin_y_for_x(x, y, bin_width)
%average of Y for each bin of X will be associated with the lower bound of the X range

ymean_binned = NaN(ceil((max(x)+abs(min(x)))/bin_width), 1);
yer_binned = NaN(ceil((max(x)+abs(min(x)))/bin_width), 1);
x_binned = NaN(ceil((max(x)+abs(min(x)))/bin_width), 1);

full_range = linspace(min(x),max(x),ceil((max(x)+abs(min(x)))/bin_width));

for b = 2:ceil((max(x)+abs(min(x)))/bin_width)
   
   range = [full_range(b-1) full_range(b)];
   
   index = find(x > range(1) & x <= range(2));
      
   ymean = mean(y(index), 'omitnan');
   yer = mean(y(index), 'omitnan')./sqrt(sum(~isnan(y(index)))); 
   
   ymean_binned(b, 1) = ymean;
   yer_binned(b, 1) = yer;
   x_bin = range(2);
   
   x_binned(b, 1) = x_bin;
   
end

end