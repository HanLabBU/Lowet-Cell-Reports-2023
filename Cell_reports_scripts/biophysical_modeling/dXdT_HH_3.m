function [f,varargout] = dXdT_HH3(~,x,I_app,I_app2, modB)
% FUNCTION dXdT_HH
% Inputs: t - time (milliseconds)
% x - vector of state variables {v,m,n,h}
% I_app - applied current (microA cm^{-2})
%
% Outputs: f - vector of time derivatives
% {dv/dt,dm/dt,dn/dt,dh/dt}
% Resting potentials, conductivities, and capacitance:
c1=1;
c2=c1;
p=0.15;
V_Na = 55;
V_K = -90;
V_L = -65;

g_Na = 55;
g_K = 20;
g_L = 0.18;  %0.3
g_L2=g_L;
g_NaP=0.04+modB./10;
g_KS=0.4+modB./1;

C_m = 1;%1e-6;
C_m2=C_m;
% State Variables:
v = x(1);
m = x(2);
n = x(3);
h = x(4);
v2 = x(5);
mm = x(6);
kk = x(7);
% alphas and betas:
a_m = -0.1*(v+31)/(exp(-0.1*(v+31))-1);
b_m = 4*exp(-(v+56)/18);
a_h = 0.07*exp(-(v+47)/20);
b_h = 1./(exp(-0.1*(v+17)) + 1);
a_n = -0.01*(v+34)./(exp(-0.1*(v+34))-1);
b_n = 0.125*exp(-(v+44)/80);
% a_m = 0.1*(25-v)/(exp((25-v)/10)-1);
% b_m = 4*exp(-v/18);
% a_h = 0.07*exp(-v/20);
% b_h = 1 ./ (exp((30-v)/10) + 1);
% a_n = 0.01*(10-v)./(exp((10-v)/10)-1);
% b_n = 0.125*exp(-v/80);
a_mm= 1/((exp(-(v2+57.7)./7.7)+1));
a_kk= 1/((1+exp(-(v2+35)./6.5)));
a_t=(200/(exp(-(v2+55)/30)+exp((v2+55)/30)));
% Computing currents:
I_Na = (m^3)*h*g_Na*(v-V_Na);
I_K = (n^4)*g_K*(v-V_K);
I_L = g_L*(v-V_L);

I_L2 = g_L2*(v2-V_L);
I_NaP = (mm^3)*g_NaP*(v2-V_Na);
I_KS = (kk)*g_KS*(v2-V_K);
% Computing derivatives:
f(1) = (-I_Na - I_K - I_L - (c1./p)*((v-v2)) + I_app )/C_m;
f(2) = 10.*(a_m*(1-m) - b_m*m);
f(3) = 3.3.*(a_n*(1-n) - b_n*n);
f(4) = 3.3.*(a_h*(1-h) - b_h*h);
f(5) = (-I_NaP - I_KS - I_L2  - (c2./(1-p))*((v2-v)) + I_app2)/C_m2;
f(6)= (1*(a_mm-mm));
f(7)= (1*(a_kk-kk))/a_t;
f=f';
%f(8)= a_kk*(1-kk);
% Outputting the conductivities
varargout{1} = [(m^3)*h*g_Na (n^4)*g_K g_L];