clear all;
clc;

%% Constants
global b h M C_s g C_w C_p Ix Iy Iz
b = 0.15; %m
h = 0.05; % m
M = 1.25; %kg
C_s = 0.85;
g =9.8;
C_w = 0.25;
C_p = 0.9;


% Moment of Inertia ellipsoid
ae = b; be = b/1.5; ce = be;
Iz = 0.2*M*(ae^2 + be^2);
Ix = 0.2*M*(ce^2 + be^2) + M*h^2;
Iy = 0.2*M*(ce^2 + ae^2) + M*h^2;
const = [Ix;Iy;Iz;h;b;M;g;C_s;C_p;C_w];

%% Dyanmics -------------------------------

% dx/dt = f(x,u)
% x = [u, theta, d_theta, psi, d_psi]
% U = Tau;     [A, Omega]
n = 4; % no. of states
m = 1;  % no. of inputs
f_xu = @eom_grnd_roll_CS_red; % f(x,u) -- non-linear dynamics with 1 input torque

%% Discretize -----------------------------
dt = 0.001;

% Runge-Kutta 4
k1 = @(t,x,u) (f_xu(t,x,const,u));
k2 = @(t,x,u) (f_xu(t + dt/2,x + k1(t,x,u)*dt/2,const,u));
k3 = @(t,x,u) (f_xu(t + dt/2,x + k2(t,x,u)*dt/2,const,u));
k4 = @(t,x,u) (f_xu(t + dt,x + k3(t,x,u)*dt,const,u));
f_ud = @(t,x,u) (x + (dt/6)*( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)));

%% Data Collection ------------------------
Nsim = 2000;
Ntraj = 1000;

% Random forcing
Tau = 150*rand([Nsim Ntraj])-75;  % A = 80*rand([Nsim Ntraj]); % Omega = 15*rand([Nsim Ntraj]);

% Random initial conditions
X01 = (rand(1,Ntraj)*25); X03 = 2*pi*(rand(1,Ntraj)*2-1);
X02 = 40*(rand(1,Ntraj)*2-1); X04 = (rand([1,Ntraj])*2-1)*20;
Xcurrent = [X01;X02(1,:);X03(1,:);X04(1,:)];

X = []; Y = []; U = [];
for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,Tau(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Tau(i,:)];
    Xcurrent = Xnext;
end
%% Basis (lifting) functions --------------
% cubic monomial for 4 state variable has 56 basis functions
tic
% vn = size(X,1);
% syms y [vn 1]
syms u theta_dot psi psi_dot
x = [u;theta_dot;psi;psi_dot];
zm = monomials(x,[0:4]);
parfor i = 1:size(X,2)
    i
    Xlift(:,i) = double(subs(zm,x,X(:,i)));
    Ylift(:,i) = double(subs(zm,x,Y(:,i)));
    clc;
end
et_data_lift = toc;

%% Linear predictor -----------------------
Nlift = size(zm,1);
W = [Ylift ; X];
V = [Xlift; U];
VVt = V*V';
WVt = W*V';
Mg = WVt * pinv(VVt); % Matrix [A B; C 0]
Alift = Mg(1:Nlift,1:Nlift);
Blift = Mg(1:Nlift,Nlift+1:end);
Clift = Mg(Nlift+1:end,1:Nlift);
save('koopman_roll_CS_not_unif.mat')