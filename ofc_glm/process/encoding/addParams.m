function [wMLnew,opsnew] =  addParams(wMLold,params,index,ops,sizeops)
    %will insert params and update the ops struct
    
    %sizing info that will be added to ops
    ls = sizeops.ls; %kernel length
    dimst = sizeops.dimst;
    dimbcst = sizeops.dimbcst;
    nkern = ops.nkern+1;
    
    
    opsnew = ops;
    opsnew.nkern = nkern;
    opsnew.ls = [ops.ls(1:index-1),ls,ops.ls(index:end)];
    opsnew.dimst = [ops.dimst(1:index-1),dimst,ops.dimst(index:end)];
    opsnew.dimbcst = [ops.dimbcst(1:index-1),dimbcst,ops.dimbcst(index:end)];
    
    paramstart = sum(ops.dimst(1:index-1));
    
    wMLnew = [wMLold(1:paramstart),params,wMLold(paramstart+1:end)];
   
    
    
    