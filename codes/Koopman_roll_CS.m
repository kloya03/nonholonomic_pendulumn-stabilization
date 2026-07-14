clear all;
clc;
tic

%% Constants
global b h M C_s g C_w C_p Ix Iy Iz
b = 0.15; %m 
h = 0.05; % m
M = 2.5; %kg
C_s = 0.85;
g =9.8;
C_w = 0.85;
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
n = 5; % no. of states
m = 1;  % no. of inputs
f_xu = @eom_grnd_roll_CS; % f(x,u) -- non-linear dynamics with 1 input torque

%% Discretize -----------------------------
dt = 0.01;

% Runge-Kutta 4
k1 = @(t,x,u) (f_xu(t,x,const,u));
k2 = @(t,x,u) (f_xu(t + dt/2,x + k1(t,x,u)*dt/2,const,u));
k3 = @(t,x,u) (f_xu(t + dt/2,x + k2(t,x,u)*dt/2,const,u));
k4 = @(t,x,u) (f_xu(t + dt,x + k3(t,x,u)*dt,const,u));
f_ud = @(t,x,u) (x + (dt/6)*( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)));

%% Data Collection ------------------------
Nsim = 200;
Ntraj = 1000;

% Random forcing
Tau = 150*rand([Nsim Ntraj])-75;  % A = 80*rand([Nsim Ntraj]); % Omega = 15*rand([Nsim Ntraj]);

% Random initial conditions
X01 = (rand(1,Ntraj)*25); X02 = 2*pi*(rand(2,Ntraj)*2-1);
X03 = 20*(rand(2,Ntraj)*2-1);
Xcurrent = [X01;X02(1,:);X03(1,:);X02(2,:);X03(2,:)];

X = []; Y = []; U = [];
for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,Tau(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Tau(i,:)];
    Xcurrent = Xnext;
end

%% Basis (lifting) functions --------------
% cubic monomial for 5 state variable has 56 basis functions
tic
vn = size(X,1);
syms y [vn 1]
zm = monomials(y,[0:3]);
parfor i = 1:size(X,2)
    Xlift(:,i) = double(subs(zm,y,X(:,i)));
    Ylift(:,i) = double(subs(zm,y,Y(:,i)));
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
save('koopman_roll_CS.mat')

% %% Simulation Comparison -----------------
% Tmax = 10;
% nsim = Tmax/dt;
% u_dt = @(i)(-75*sin((i*dt))); % control signal
% 
% % Initial condition
% x0 = [0.1;0.1;0.02;0.01;0.01];
% x_true = x0;
% 
% % Lifted initial condition
% 
% xlift = double(subs(zm,y,x0));
% 
% % Simulate
% for i = 0:nsim-1
%     x_liftt = double(subs(zm,y,x_true(:,end)));
%     % Koopman predictor
% %     xlift = [xlift, Alift*xlift(:,end) + Blift*u_dt(i)]; % Lifted dynamics
%     xlift = [xlift, Alift*x_liftt + Blift*u_dt(i)]; % Lifted dynamics
% 
%     % True dynamics
%     x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
% end
% x_koop = Clift * xlift; % Koopman predictions
% et = toc;
% % [t1,sts1] = ode45(@(t,sts)eom_grnd_roll_CS_input(t,sts,const,inp),[0 10],x0);
% 
% %% Plotting -----------------------------
% str1 = sprintf('Koopman simulation with %1.0f sim*%1.0f traj',Nsim,Ntraj);
% lw=2;
% figure(1)
% nexttile
% plot([0:nsim-1]*dt,u_dt(0:nsim-1),'linewidth',lw); hold on
% ylabel('torque input ', 'interpreter','latex');
% xlabel('Time [s]','interpreter','latex');
% set(gca,'fontsize',20)
% 
% nexttile
% plot([0:nsim]*dt,x_true(1,:),'linewidth',lw); hold on
% plot([0:nsim]*dt,x_koop(1,:), '--r','linewidth',lw)
% % axis([0 Tmax min(x_koop(1,:))-0.15 max(x_koop(1,:))+0.15])
% title(str1)
% ylabel('$u$','interpreter','latex');
% xlabel('Time [s]','interpreter','latex');
% set(gca,'fontsize',20)
% LEG = legend('True RK-4','Koopman','ODE45','location','south');
% set(LEG,'interpreter','latex')
% 
% nexttile
% plot([0:nsim]*dt,x_true(2,:),'linewidth',lw); hold on
% plot([0:nsim]*dt,x_koop(2,:), '--r','linewidth',lw)
% % axis([0 Tmax min(x_koop(2,:))-0.15 max(x_koop(2,:))+0.15])
% ylabel('$\theta$ yaw','interpreter','latex');
% xlabel('Time [s]','interpreter','latex');
% set(gca,'fontsize',20)
% % LEG = legend('True','Koopman','location','southwest');
% % set(LEG,'interpreter','latex')
% 
% nexttile
% plot([0:nsim]*dt,x_true(3,:),'linewidth',lw); hold on
% plot([0:nsim]*dt,x_koop(3,:), '--r','linewidth',lw)
% % axis([0 Tmax min(x_koop(3,:))-0.15 max(x_koop(3,:))+0.15])
% ylabel('$\dot{\theta}$ yaw rate','interpreter','latex');
% xlabel('Time [s]','interpreter','latex');
% set(gca,'fontsize',20)
% % LEG = legend('True','Koopman','location','southwest');
% % set(LEG,'interpreter','latex')
% 
% nexttile
% plot([0:nsim]*dt,x_true(4,:),'linewidth',lw); hold on
% plot([0:nsim]*dt,x_koop(4,:), '--r','linewidth',lw)
% % axis([0 Tmax min(x_koop(4,:))-0.15 max(x_koop(4,:))+0.15])
% ylabel('$\psi$ roll','interpreter','latex');
% xlabel('Time [s]','interpreter','latex');
% set(gca,'fontsize',20)
% % LEG = legend('True','Koopman','location','southwest');
% % set(LEG,'interpreter','latex')
% 
% nexttile
% plot([0:nsim]*dt,x_true(5,:),'linewidth',lw); hold on
% plot([0:nsim]*dt,x_koop(5,:), '--r','linewidth',lw)
% % axis([0 Tmax min(x_koop(5,:))-0.15 max(x_koop(5,:))+0.15])
% ylabel('$\dot{\psi}$ roll rate','interpreter','latex');
% xlabel('Time [s]','interpreter','latex');
% set(gca,'fontsize',20)
% % LEG = legend('True','Koopman','location','southwest');
% % set(LEG,'interpreter','latex')
% saveas(figure(1),'Predictor_comparison')
