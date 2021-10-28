function [basefun] = chooseBasis(type,tmesh,params)
%CHOOSEBASIS: create basis functions for spike and stimulus filters
%
%
%input: 
%   %type: a handle to decide which basis to use
    %tmesh: vector of time points. just needs to be 1:numel(tmesh)
    %params: the importnat parameters for each model. if not
%            supplied, defaults found by inspection given below
%            form: params = [lh,dim,dimbc,ls,dimst,dimbcst];
%                   lh: length history filter
%                   dim: number of basis fun. coponents of history filter
%                   dimbc: number of boxcar funs as part of history filter
%                   ls: length stimulus filter
%                   dimst: number of basis fun. coponents of stim filter
%                   dimbcst: number of boxcar funs as part of stim filter

%outputs:
%basefun: [dim x ls] matrix of basis functions for history filter

%update: 20200811. adding extra form (temporarily) in which boxcar
%functions are decyaing exponentials over a few time bins

Nt = numel(tmesh);

%new basis, 100 ms each
switch type
    case 'none'
        %each timepoint is a parameter
        ls = params(1); %length history filter
        dim = params(2); %number of history basis func
        dimbc = params(3); %number of boxcar funcion, history

        %phivecst = phivec;
        basefun = eye(ls);


    case 'gaussian'
        %tile the space with gaussians
        ls = params(1); %length history filter
        dim = params(2); %number of history basis func
        dimbc = params(3); %number of boxcar funcion, history

        %set width of gaussian to overlap at 0.5
        sigma_h = ls/(2*dim*sqrt(-log(0.5)));

        %find locations of each gaussian
        basefun_loc = linspace(1,ls,dim);
        basefun = zeros(dim,ls);

        lsvec = 1:ls;
        for j = 1:dim
            basefun(j,:) = exp(-(lsvec-basefun_loc(j)).^2/sigma_h^2);
        end
        
    case 'lip'
        %LIP-type basis functions: with boxcars for short times
        %and raised log-scaled cosine for others
        %params = [lh,dim,dimbc,ls,dimst,dimbcst];
        if isempty(params) || nargin <= 2
            params = [199,11,5,100,11,5];       %old params    
        else
            ls = params(1); %length history filter
            dim = params(2); %number of history basis func
            dimbc = params(3); %number of boxcar funcion, history
        end

        wbc = 2; %width of boxcar
        
        %or manually create basis functions with different sizes.
        %params tuned by inspection
        c = 0.9;
        a = 3.3;
        phivec = linspace(7,12,dim-dimbc);
        basefun = weber_basis(tmesh,a,c,dim-dimbc,phivec);
        
        %trim edges of basefun to not interfere with boxcars
        basefun(:,1:dimbc*wbc)=0;

        %add boxcar functions to beginning of filter. used to
        %capture short refractory period features
        basebc = zeros(dimbc,Nt);
        for j = 1:dimbc
            basebc(j,(j-1)*wbc+1:j*wbc) = 1;
        end

        basefun = [basebc;basefun];
        %trim to correct size
        basefun = basefun(:,1:ls);

    case 'ofc'
        %basis functions for ofc neurons: with boxcars for short times
        %and raised log-scaled cosine for others
        %params = [lh,dim,dimbc,ls,dimst,dimbcst];
        if isempty(params) || nargin <= 2
            params = [50,5,2,50,5,2];           
        else
            ls = params(1); %length kernel
            dim = params(2); %number of history basis func
            dimbc = params(3); %number of boxcar functions
        end

        wbc = 1; %width of boxcar
        
        % manually create basis functions with different sizes.
        %params tuned by inspection        
        c = 0.9;
        a = 3.8; %good for lh = 50;
        ast = 2.3;
        
        phivec = linspace(7,12,dim-dimbc);
        basefun = weber_basis(tmesh,a,c,dim-dimbc,phivec);

        %trim edges of basefun to not interfere with boxcars
        basefun(:,1:dimbc*wbc)=0;

        %add boxcar functions to beginning of filter. used to
        %capture short refractory period features
        basebc = zeros(dimbc,Nt);
        for j = 1:dimbc
            basebc(j,(j-1)*wbc+1:j*wbc) = 1;
        end
        basefun = [basebc;basefun];

        %trim to correct size
        basefun = basefun(:,1:ls);

    case 'ofc_adaptive'
        %to try and chnage the coeff a and phi vector based on
        %number of basis functions and filter length

        ls = params(1); %length kernel
        dim = params(2); %number of history basis func
        dimbc = params(3); %number of boxcar functions

        wbc = 1; %width of boxcar         
        c = 0; %shift parameter. not used 
        
        %check for if dimensions mean turn off one params set
        if dim > 0

            %evenly tile the functions from beginning of filt to ~45% of history filter
            xrange = floor([1,0.45*ls]); %min and max crossing points
            crossvec = floor(linspace(xrange(1),xrange(2),dim-dimbc+1)); %locations of mid points
            maxvec = floor(diff(crossvec)/2+crossvec(1:end-1)); %locations of maxima

            %set scaling by describing smoothness of first function
            %idea: cos(a*log(tmax_1)-phi_1) = 1. Find what angle of rotation is
            %needed to set a nearby point (crossvec(2)=halway between first and
            %second maxima) to a given value. implies
            %a*log(tmax_1)-a*log(crossvec(2)) = theta_cross. 
            %backsolve for a

            crossy = 0.2; %height at which first basis fun should hit on right side
            theta_cross = acos(2*crossy-1); %relative angle change required from maxima

            rcross = crossvec(2)/maxvec(1);
            a = (theta_cross)/log(rcross);       
            thetavec = a*log(maxvec);
            phivec = thetavec;
            
            basefun = weber_basis(tmesh,a,c,dim-dimbc,phivec);
            
        else
            basefun = [];
        end      

        %trim edges of basefun to not interfere with boxcars
        %offset from first element for spike history, as it gets removed
        %from filter() and would lead to rank-deficiency of Hessian of -LL
        basefun(:,2:dimbc*wbc+1)=0; 

        %add boxcar functions to beginning of filter. used to
        %capture short refractory period features
        basebc = zeros(dimbc,Nt);
        for j = 1:dimbc
            basebc(j,(j-1)*wbc+2:j*wbc+1) = 1;
        end
        basefun = [basebc;basefun];

        %trim to correct size
        basefun = basefun(:,1:ls);
        
    case 'ofc_adaptive_logUniform'
        %basis functions are uniformly spaced in log time. so more at
        %beginning of time, and less as time progresses.
        %basic placement of extrema are given by t_j = log(mj+b)t_Delta
        %will also place things relative to a particular window of support
        %this will be the final param
        

        ls = params(1); %length kernel
        dim = params(2); %number of history basis func
        dimbc = params(3); %number of boxcar functions
        ls_supp = params(4); %reference length of the support for placing the basis functions. 
        dim_supp = params(5); %refrence number of basis funs 
        %in general ls <= ls_supp to allow for 
        M = dim_supp-dimbc;

        wbc = 5; %width of boxcar. 5 for exp decay. 1 for impulse       
        c = 0; %shift parameter. not used 
                
       
        
        %check for if dimensions mean turn off one params set
        if dim_supp >0
       
            tau_s = floor(0.40*ls_supp); %final baiss function max
            tau_delta = tau_s/(M+1); %a spacing parameter
            b = tau_delta/2;
            
            %m = tau_delta*2; %2/log(2) keeps the log version sublinear 
            %maxvec = ceil(log(1:M)*m+b) %locations of maxima
            
            %try sublinear exponential scaling
            m = log(tau_delta*(M+1/2))/M;
            b = tau_delta/2-1; %-1 do account for exponential term ~1
            maxvec = ceil(exp((1:M)*m)+b);

            %set scaling by describing smoothness of first function
            %idea: cos(a*log(tmax_1)-phi_1) = 1. Find what angle of rotation is
            %needed to set a nearby point (crossvec(2)=halway between first and
            %second maxima) to a given value. implies
            %a*log(tmax_1)-a*log(crossvec(2)) = theta_cross. 
            %backsolve for a

            crossy = 0.2; %height at which first basis fun should hit on right side
            theta_cross = acos(2*crossy-1); %relative angle change required from maxima

            rcross = 2; %rcross = t_a/t_1
            a = (theta_cross)/log(rcross);       
            thetavec = a*log(maxvec);
            phivec = thetavec;
            
            basefun = weber_basis(tmesh,a,c,dim_supp-dimbc,phivec);
            
        else
            basefun = [];
        end      

        %trim edges of basefun to not interfere with boxcars
        %offset from first element for spike history, as it gets removed
        %from filter() and would lead to rank-deficiency of Hessian of -LL
        
        %basefun(:,2:dimbc*wbc+1)=0; 

        %add boxcar functions to beginning of filter. used to
        %capture short refractory period features
        basebc = zeros(dimbc,Nt);
        for j = 1:dimbc
            %basebc(j,(j-1)*wbc+1:j*wbc) = 1; %CORRECT
            %basebc(j,(j-1)*wbc+2:j*wbc+1) = 1; %wrong. keep there temporarily
            basebc(j,j:j+wbc) = exp(-(0:wbc)); %transferring decaying exponential form over
        end
        basefun = [basebc;basefun];

        %trim to correct size
        basefun = basefun(1:dim,1:ls);
        
    %TODO: REMOVE    
    case 'ofc_adaptive_logUniform_decExp'
        %as above, but with decaying exponential functions rather that
        %delta functions
        

        ls = params(1); %length kernel
        dim = params(2); %number of history basis func
        dimbc = params(3); %number of boxcar functions
        ls_supp = params(4); %reference length of the support for placing the basis functions. 
        dim_supp = params(5); %refrence number of basis funs 
        %in general ls <= ls_supp to allow for 
        M = dim_supp-dimbc;

        wbc = 5; %width of boxcar. changed        
        c = 0; %shift parameter. not used 
                
       
        
        %check for if dimensions mean turn off one params set
        if dim_supp >0
       
            tau_s = floor(0.40*ls_supp); %final baiss function max
            tau_delta = tau_s/(M+1); %a spacing parameter
            b = tau_delta/2;
            
            %m = tau_delta*2; %2/log(2) keeps the log version sublinear 
            %maxvec = ceil(log(1:M)*m+b) %locations of maxima
            
            %try sublinear exponential scaling
            m = log(tau_delta*(M+1/2))/M;
            b = tau_delta/2-1; %-1 do account for exponential term ~1
            maxvec = ceil(exp((1:M)*m)+b);

            %set scaling by describing smoothness of first function
            %idea: cos(a*log(tmax_1)-phi_1) = 1. Find what angle of rotation is
            %needed to set a nearby point (crossvec(2)=halway between first and
            %second maxima) to a given value. implies
            %a*log(tmax_1)-a*log(crossvec(2)) = theta_cross. 
            %backsolve for a

            crossy = 0.2; %height at which first basis fun should hit on right side
            theta_cross = acos(2*crossy-1); %relative angle change required from maxima

            rcross = 2; %rcross = t_a/t_1
            a = (theta_cross)/log(rcross);       
            thetavec = a*log(maxvec);
            phivec = thetavec;
            
            basefun = weber_basis(tmesh,a,c,dim_supp-dimbc,phivec);
            
        else
            basefun = [];
        end      

        %trim edges of basefun to not interfere with boxcars
        %offset from first element for spike history, as it gets removed
        %from filter() and would lead to rank-deficiency of Hessian of -LL
        %basefun(:,2:dimbc*wbc+1)=0; 

        %add boxcar functions to beginning of filter. used to
        %capture short refractory period features
        basebc = zeros(dimbc,Nt);
        for j = 1:dimbc
            %basebc(j,(j-1)*wbc+1:j*wbc) = 1; original delta fun
            basebc(j,j:j+wbc) = exp(-(0:wbc));
        end
        basefun = [basebc;basefun];

        %trim to correct size
        basefun = basefun(1:dim,1:ls);
        
        
    %TODO: remove    
    case 'ofc_adaptive_plus'
        %same as ofc_adaptive_loguniform, with just an extra few at a few
        %seconds later. THIS IS A CUSTOM BASIS
        

        ls = params(1); %length kernel
        dim = params(2); %number of history basis func
        dimbc = params(3); %number of boxcar functions
        ls_supp = params(4); %reference length of the support for placing the basis functions. 
        dim_supp = params(5); %refrence number of basis funs 
        %in general ls <= ls_supp to allow for 
        M = dim_supp-dimbc;

        wbc = 5; %width of boxcar. 5 for exp decay. 1 for impulse       
        c = 0; %shift parameter. not used 
                
       
        
        %check for if dimensions mean turn off one params set
        if dim_supp >0
       
            tau_s = floor(0.40*ls_supp); %final baiss function max
            tau_delta = tau_s/(M+1); %a spacing parameter
            b = tau_delta/2;
            
            %m = tau_delta*2; %2/log(2) keeps the log version sublinear 
            %maxvec = ceil(log(1:M)*m+b) %locations of maxima
            
            %try sublinear exponential scaling
            m = log(tau_delta*(M+1/2))/M;
            b = tau_delta/2-1; %-1 do account for exponential term ~1
            maxvec = ceil(exp((1:M)*m)+b);

            %set scaling by describing smoothness of first function
            %idea: cos(a*log(tmax_1)-phi_1) = 1. Find what angle of rotation is
            %needed to set a nearby point (crossvec(2)=halway between first and
            %second maxima) to a given value. implies
            %a*log(tmax_1)-a*log(crossvec(2)) = theta_cross. 
            %backsolve for a

            crossy = 0.2; %height at which first basis fun should hit on right side
            theta_cross = acos(2*crossy-1); %relative angle change required from maxima

            rcross = 2; %rcross = t_a/t_1
            a = (theta_cross)/log(rcross);       
            thetavec = a*log(maxvec);
            phivec = thetavec;
            
            basefun = weber_basis(tmesh,a,c,dim_supp-dimbc,phivec);
            
        else
            basefun = [];
        end      

        %trim edges of basefun to not interfere with boxcars
        %offset from first element for spike history, as it gets removed
        %from filter() and would lead to rank-deficiency of Hessian of -LL
        
        %basefun(:,2:dimbc*wbc+1)=0; 

        %add boxcar functions to beginning of filter. used to
        %capture short refractory period features
        basebc = zeros(dimbc,Nt);
        for j = 1:dimbc
            %basebc(j,(j-1)*wbc+1:j*wbc) = 1; %CORRECT
            %basebc(j,(j-1)*wbc+2:j*wbc+1) = 1; %wrong. keep there temporarily
            basebc(j,j:j+wbc) = exp(-(0:wbc)); %transferring decaying exponential form over
        end
        basefun = [basebc;basefun];

        %trim to correct size
        basefun_temp = basefun(1:dim,1:ls);
        
        %add dims 5:7 over by 30 points
        nextra = 3;
        if size(basefun_temp,1) > 0
            basefun_extra = zeros(3,ls);
            for j = 1:nextra
                ind = 40+(j-1)*10;
                basefun_extra(j,ind+1:ind+20) = basefun_temp(6,4:23);
            end

            basefun = [basefun_temp;basefun_extra];
        else %no filters (typically for history filter_
            basefun = basefun_temp;
        end

        
    case 'power' %TODO: was this ever debugged?
        %power series
        ls = params(1); %length kernel
        dim = params(2); %number of terms                     
        
        fx = @(x,d) x.^(d)'/(x(end)^d); %give some normalization
        basefun = zeros(dim,ls);
        for j = 1:dim
            basefun(j,:) = fx(1:ls,j);
        end

        
        
end