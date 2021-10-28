function [basfun,a,c] = weber_basis(tmesh,a,c,dim,phivec,varargin)
%WEBER_BASIS: log-scaled cosine basis function form. weber arxiv 2016 paper 
%
%creates basis functions of form m_j(t) =  cos( a log(t+c)-phi_j)  )
%
%inputs
%tmesh: time vector
%a:     basis fun param
%c:     basis fun param
%dim:   number of basis functions
%phivec: vector of phaseshift terms for each basis function
%varargin:  timescale of basis functions

%output
%basefun: [dim x lt] matrix of basis functions, where lt = length(tmesh)
%if no params a and c supplied, but a general timescale given, then
%a: calculated amplitude parameter
%c: calculated phase parameter

if isempty(phivec)
    phivec = linspace(0,pi,dim);
end
lt = length(tmesh);
basfun = zeros(dim,lt);


% alternative mode: fit the basis function to a given timescale
if ~isempty(varargin)>0
    tscale = varargin{1}/2;
    c = 1-tscale;
    a = -pi/log(1-tscale);
    disp('a,c');
    disp([a,c])
end


for j = 1:dim
    if  length(a) > 1
        basfun(j,:) = a(j)*log(tmesh+c(j))-phivec(j);
    else
        basfun(j,:) = a*log(tmesh+c)-phivec(j);
    end
    %clear out anything outside the range [-pi,pi]
    basfun(j,basfun(j,:) < -pi | basfun(j,:) > pi) = -pi; 
    basfun(j,:) = 0.5*cos(basfun(j,:))+0.5;

end


