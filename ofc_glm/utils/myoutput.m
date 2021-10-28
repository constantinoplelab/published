function stop = myoutput(x,optimvalues,state,ops)
%   MYOUTPUT: a function to call during optimization for periodic saves
%
%   ops must be a struct containing "savefreq", number of iterations
%   between saves
% 
    savefreq = ops.savefreq;
    stop = false;
    iteration = optimvalues.iteration;
    if rem(iteration,savefreq)==0 & iteration > 0
        save(strcat([ops.fname,'.tmp']),'x','iteration')
    end
