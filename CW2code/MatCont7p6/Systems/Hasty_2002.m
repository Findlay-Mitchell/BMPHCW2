function out = Hasty_2002
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
function dydt = fun_eval(t,kmrgd,par_s,par_a,par_f,par_g,par_T)
dydt=[(1+kmrgd(1)^2+par_a*par_s*kmrgd(1)^4)/((1+kmrgd(1)^2+par_s*kmrgd(1)^4)*(1+kmrgd(2)^4))-par_f*kmrgd(1);
((1+kmrgd(1)^2+par_a*par_s*kmrgd(1)^4)/((1+kmrgd(1)^2+par_s*kmrgd(1)^4)*(1+kmrgd(2)^4))-par_g*kmrgd(2))/par_T;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Hasty_2002);
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_s,par_a,par_f,par_g,par_T)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_s,par_a,par_f,par_g,par_T)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_s,par_a,par_f,par_g,par_T)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_s,par_a,par_f,par_g,par_T)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_s,par_a,par_f,par_g,par_T)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_s,par_a,par_f,par_g,par_T)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_s,par_a,par_f,par_g,par_T)
