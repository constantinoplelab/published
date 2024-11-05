function effectsize = effsize_ind(x, y) 

meandiff = mean(x, 'omitnan') - mean(y, 'omitnan');
std_pooled = std(x, 'omitnan'); %independent groups
effectsize = meandiff./std_pooled;

end