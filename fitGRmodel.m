function [grfun, domain, gof] = fitGRmodel(t, od, varargin)
% UIFITGRMODEL fits a piecewise linear function to the exponential
% growth rate of the growth curve data in T and OD.
%
% 20141013
parser = inputParser;
addParamValue(parser,'keepplot',0,@islogical);
addParamValue(parser,'prevfit',[],@(x) isa(x,'cfit'));

parse(parser,varargin{:});
keepplot = parser.Results.keepplot;
prevfit = parser.Results.prevfit;

if isempty(prevfit)
    % get parameter guesses from user
    dod = diff(od)./diff(t');
    [hf ha] = gridplot(1,1,600,450);
    plot(t(2:end), dod, '-','color',[.8 .8 .8]);
    hold all
    
    dodsm = smooth(dod,10);
    plot(t(2:end),dodsm,'-','color',[.5 .5 .5]);
    
    xlim([0 20]);
    ylim([-.5 1]);
    plot(xlim,[0 0],'k:');
    
    disp('Click on a guess for (tau1, r1)...');
    [x,y]=ginput(1);
    r1 = y;
    tau1 = x;
    plot([0 x],[r1 r1],'k--')
    plot([x],[r1],'ko')
    
    disp('Click on a guess for (tau1 + tau2, r_min)...');
    [x,y]=ginput(1);
    r_min = y;
    tau2 = x - tau1;
    plot([tau1 tau1+tau2], [r1 r_min],'k--')
    plot([tau1+tau2], [r_min],'ko')
    
    disp('Click on a guess for (tau1 + tau2 + tau3, r2)...');
    [x,y]=ginput(1);
    r2 = y;
    tau3 = x - tau1 - tau2;
    plot([tau1+tau2 tau1+tau2+tau3], [r_min r2],'k--')
    plot([tau1+tau2+tau3], [r2],'ko')
    
    disp('Click on a guess for tau1 + tau2 + tau3 + tau4 ...');
    [x,y]=ginput(1);
    tau4 = x - tau1 - tau2 - tau3;
    plot([tau1+tau2+tau3 x], [r2 r2],'k--')
    plot([x], [r2],'ko')
    
    y0 = -7;    % assume this is where od is aligned
else
    % use parameters from a previous fit
    r1 = prevfit.r1;
    if any(strcmp(coeffnames(prevfit),'tau3'))
        r_min = prevfit.r_min;
        tau3 = prevfit.tau3;
    else
        r_min = 0;
        tau3 = 0;
    end
    r2 = prevfit.r2;
    tau1 = prevfit.tau1;
    tau2 = prevfit.tau2;
    tau4 = prevfit.tau4;
    y0 = prevfit.y0;
end

% use smoothing spline to find saturation point
odfun = csaps(t,od,0.7);
dodfun = fnder(odfun,1);
z = fnzeros(dodfun,[10 20]);
t_end = z(end);
domain = [0 t_end];

if isempty(prevfit)
    plot([tau1+tau2+tau3+tau4 t_end], [r2 0],'k--')
    plot([t_end], [0],'ko')
    
    if ~keepplot
        close(hf);
    end
end

% 4-phase fit
fitfunc = fittype(@(r1, r_min, r2, tau1, tau2, tau3, tau4, y0, x) ...
    GRmodel(r1, r_min, r2, tau1, tau2, tau3, tau4, t_end, y0, x));

idx = t >= domain(1) & t <= domain(2);
[grfun, gof] = fit(t(idx)',od(idx),fitfunc,...
    'startpoint',[r1 r_min r2 tau1 tau2 tau3 tau4 y0],...
    'lower',[0 -Inf 0 0 0 0 0 -Inf],...
    'robust','bisquare');
% gof

% 3-phase fit
fitfunc = fittype(@(r1, r2, tau1, tau2, tau4, y0, x) ...
    GRmodel(r1, [], r2, tau1, tau2, 0, tau4, t_end, y0, x));

idx = t >= domain(1) & t <= domain(2);
[grfun2, gof2] = fit(t(idx)',od(idx),fitfunc,...
    'startpoint',[r1 r2 tau1 tau2 tau4 y0],...
    'lower',[0 0 0 0 0 -Inf],...
    'robust','bisquare');
% gof2

% return the better fit
if gof2.sse < gof.sse
    grfun = grfun2;
    gof = gof2;
end

end

function y = GRmodel(r1, r_min, r2, tau1, tau2, tau3, tau4, t_end, y0, t)
% computes a diauxic growth curve using geometrically defined parameters
%
% 20141014

y = nan(length(t),1);

t1 = tau1;
idx = t >= 0 & t < t1;
y(idx) = y0 + r1.*t(idx);
y1 = y0 + r1.*t1;

if tau3 > 0
    % 4-phase model
    t2 = t1 + tau2;
    idx = t >= t1 & t < t2;
    int1 = @(t) (r1.*t + 1/2.*(r_min-r1)./tau2.*t.^2 ...
        - (r_min-r1).*t1./tau2.*t);
    y(idx) = y1 + int1(t(idx)) - int1(t1);
    y2 = y1 + int1(t2) - int1(t1);
    
    t3 = t2 + tau3;
    idx = t >= t2 & t < t3;
    int2 = @(t) (r_min.*t + 1/2.*(r2 - r_min)./tau3.*t.^2 ...
        - (r2 - r_min).*t2./tau3.*t);
    y(idx) = y2 + int2(t(idx)) - int2(t2);
    y3 = y2 + int2(t3) - int2(t2);
    
    t4 = t3 + tau4;
    idx = t >= t3 & t < t4;
    y(idx) = y3 + r2.*(t(idx) - t3);
    y4 = y3 + r2.*(t4 - t3);
else
    % 3-phase model
    t2 = tau1 + tau2;
    idx = t >= t1 & t < t2;
    int1 = @(t) (r1.*t + 1/2.*(r2-r1)./tau2.*t.^2 ...
        - (r2-r1).*t1./tau2.*t);
    y(idx) = y1 + int1(t(idx)) - int1(t1);
    y2 = y1 + int1(t2) - int1(t1);
    
    t4 = tau1 + tau2 + tau4;
    idx = t >= t2 & t < t4;
    y(idx) = y2 + r2.*(t(idx) - t2);
    y4 = y2 + r2.*(t4 - t2);
end

t5 = t_end;
tau5 = t_end - t4;
idx = t >= t4 & t <= t5;
int4 = @(t) (r2.*t - 1/2.*r2./tau5.*t.^2 + r2./tau5.*t4.*t);
y(idx) = y4 + int4(t(idx)) - int4(t4);

% y(isnan(y)) = 0;
% out = odbg + od0.*out;

end