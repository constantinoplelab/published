function [] = dispL(s,Lmin,L)
%dispL: displays at different levels of verbosity
%if suppled L >= Lmin, message is displayed. otherwise not displayed
%allows for differnet levels of output

if L >= Lmin
    disp(s)
end