%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code: Lorenz
%
% To solve the chaotic Lorenz
% equations using 4th order
% Runge-Kutta
%
% (taken from https://en.wikipedia.org/wiki/Lorenz_system)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X=lorenzz(VAR,N)
if(nargin<1 || isempty(VAR)), VAR=10; end
if(nargin<2 || isempty(N)), N=40; end


sigma=VAR;
beta=8/3;
rho=28;
f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
%'f' is the set of differential equations and 'a' is an array containing values of x,y, and z variables.
%'t' is the time variable
[t,a] = ode45(f,linspace(0,20,N),[1 1 1]);%'ode45' uses adaptive Runge-Kutta method of 4th and
                                 %5th order to solve differential equations

X=a;
return;


set_figure(20); 
clf; plot3(a(:,1),a(:,2),a(:,3)) %'plot3' is the command to make 3D plot


set_figure(19); hold all;
%plot(t,y1,'r--');
plot(t,X(:,2),'k');



