function out = plotGRmodelfit(grfun,domain)
% PLOTGRMODELFIT plots the results of a model fit and returns the figure
% handle of the first segment
% 
% 20141014
t1 = grfun.tau1;
out = plot([domain(1) t1], [grfun.r1 grfun.r1],'-','color','r');

if any(strcmp(coeffnames(grfun),'tau3'))
    t2 = grfun.tau1 + grfun.tau2;
    plot([t1 t2], [grfun.r1 grfun.r_min],'-','color',[0 .75 0]);
    t3 = t2 + grfun.tau3;
    plot([t2 t3], [grfun.r_min grfun.r2],'-','color','b');
else
    t2 = grfun.tau1 + grfun.tau2;
    plot([t1 t2], [grfun.r1 grfun.r2],'-','color',[0 .75 0]);
    t3 = t2;
end

t4 = t3 + grfun.tau4;
plot([t3 t4], [grfun.r2 grfun.r2],'-','color','r');

if any(strcmp(coeffnames(grfun),'tau5'))
    plot([t4 t4+grfun.tau5], [grfun.r2 0],'-','color',[0 .75 0]);
else
    plot([t4 domain(2)], [grfun.r2 0],'-','color',[0 .75 0]);
end
plot(domain,[0 0],'k--');
end