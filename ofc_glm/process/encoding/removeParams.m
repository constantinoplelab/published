function [wMLnew,opsnew] =  removeParams(wMLold,index,ops)
    %will remove params and update the ops struct
    
    
    
    opsnew = ops;
    opsnew.nkern = ops.nkern-1;
    opsnew.ls = [ops.ls(1:index-1),ops.ls(index+1:end)];
    opsnew.dimst = [ops.dimst(1:index-1),ops.dimst(index+1:end)];
    opsnew.dimbcst = [ops.dimbcst(1:index-1),ops.dimbcst(index+1:end)];
    
    paramstart = sum(ops.dimst(1:index-1));
    paramstart2 = sum(ops.dimst(1:index));
    
    wMLnew = [wMLold(1:paramstart),wMLold(paramstart2+1:end)];
   
    
    
    