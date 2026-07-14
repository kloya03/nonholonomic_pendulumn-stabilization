clear all;
clc;

%% Constants
global b h M C_s g C_w C_p Ix Iy Iz
b = 0.15; %m
h = 0.05; % m
M = 1.5; %kg
C_s = 0.85;
g =9.8;
C_w = 0.5;
C_p = 0.9;


% Moment of Inertia ellipsoid
ae = b; be = b/1.5; ce = be;
Iz = 0.2*M*(ae^2 + be^2);
Ix = 0.2*M*(ce^2 + be^2) + M*h^2;
Iy = 0.2*M*(ce^2 + ae^2) + M*h^2;
const = [Ix;Iy;Iz;h;b;M;g;C_s;C_p;C_w];
Om = 25;
Am = 45;
syms Aw Bw Bu Au Uc

f = [Aw^2*b*M + Bw^2*b*M - 2*Uc*C_s==0;...
    4*M*Om*Bu + 2*Aw*Bw*b*M - 2*Au*C_s==0;...
    -4*M*Om*Au - Aw^2*b*M + Bw^2*b*M - 2*Bu*C_s==0;...                          function f(Xs)
    2*(Iz + M*b^2)*Om*Bw - Au*Bw*b*M + Aw*Bu*b*M - 2*Aw*Uc*b*M + 2*Am==0;...
    -2*(Iz + M*b^2)*Om*Aw - Au*Aw*b*M - Bw*Bu*b*M - 2*Bw*Uc*b*M==0];

% assume(Am,'real')
% assume(Om,'real')
% assume(Am>0)
% assume(Om>0)
S = solve(f);
double(S.Uc)
% tsim = 10;
% x0 = [0;0;0;0;0];
% inp = [double(S.Am);Om];
% [t1,sts1] = ode45(@(t,sts)eom_CS(t,sts,const,inp),[0 tsim],x0);
% 

