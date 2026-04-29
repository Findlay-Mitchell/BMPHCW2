function out = Hasty_Attempt_3
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_SIGMA,par_ALPHA,par_GAMMA_X,par_GAMMA_Y,par_TAU)
dydt=[(1+kmrgd(1)^2+par_ALPHA*par_SIGMA*kmrgd(1)^4)/((1+kmrgd(1)^2+par_SIGMA*kmrgd(1)^4)*(1+kmrgd(2)^4))-par_GAMMA_X*kmrgd(1);
((1+kmrgd(1)^2+par_ALPHA*par_SIGMA*kmrgd(1)^4)/((1+kmrgd(1)^2+par_SIGMA*kmrgd(1)^4)*(1+kmrgd(2)^4))-par_GAMMA_Y*kmrgd(2))/par_TAU;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Hasty_Attempt_3);
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_SIGMA,par_ALPHA,par_GAMMA_X,par_GAMMA_Y,par_TAU)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_SIGMA,par_ALPHA,par_GAMMA_X,par_GAMMA_Y,par_TAU)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_SIGMA,par_ALPHA,par_GAMMA_X,par_GAMMA_Y,par_TAU)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_SIGMA,par_ALPHA,par_GAMMA_X,par_GAMMA_Y,par_TAU)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_SIGMA,par_ALPHA,par_GAMMA_X,par_GAMMA_Y,par_TAU)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_SIGMA,par_ALPHA,par_GAMMA_X,par_GAMMA_Y,par_TAU)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_SIGMA,par_ALPHA,par_GAMMA_X,par_GAMMA_Y,par_TAU)
