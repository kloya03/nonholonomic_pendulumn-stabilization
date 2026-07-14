clc;
clear;
tic
addpath(genpath('/home/kloya/Documents/MATLAB/GPOPS-II-Distribution'))
% X_state = [x_pos, x_dot, y_pos, y_dot, z_pos, z_dot, phi, phi_dot, theta,...
%           theta_dot, psi, psi_dot, w1, w2, w3, w4];
% U_control = [a1, a2, a3, a4];

%%
%----------------- Constants parameters ------------------------%
%-----------------------------------------------------%
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

%%
%----------------- Variable parameters ------------------------%
%-----------------------------------------------------%
global t0 tf x0 xf uMin uMax

t0 = 0;
tf = 5;

x0 = [0 0 1.9 1.12 0.9]; % initial state
xf = [0 0 0 0 0];
uMin = -75;
uMax = 75;
%%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf;
bounds.phase.finaltime.upper = tf;
bounds.phase.initialstate.lower = x0;
bounds.phase.initialstate.upper = x0;
bounds.phase.state.lower = -100*ones(1,5);
bounds.phase.state.upper = 100*ones(1,5);
bounds.phase.finalstate.lower = -100*ones(1,5);
bounds.phase.finalstate.upper = 100*ones(1,5); % you can relax the end position
bounds.phase.control.lower = uMin;
bounds.phase.control.upper = uMax;
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 10000;
% bounds.phase.path.lower = 0;
% bounds.phase.path.upper =1e-4;
% bounds.eventgroup.lower = 0;
% bounds.eventgroup.upper = 1e-04;
%%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time = [0;1];
guess.phase.state = [x0;-0.3864 1.9000 -22.2722 0.900 21.8783];
guess.phase.control = [70;-50];
guess.phase.integral = 300;

%%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method = 'hp-PattersonRao';
mesh.tolerance = 1e-3;
mesh.colpointsmin = 30;
mesh.colpointsmax = 40;
mesh.maxiterations = 30;

%%
%------------- Assemble Information into Problem Structure ---------------%
%-------------------------------------------------------------------------%
setup.name = 'roll_stabilization_CS';
setup.scales.method = 'automatic-bounds';
setup.functions.continuous = @fContinuous ;
setup.functions.endpoint = @fEndpoint ;
setup.bounds = bounds ;
setup.guess = guess ;
setup.mesh = mesh ;
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.tolerance=1e-6; %Default is 1e-7
setup.nlp.ipoptoptions.maxiterations = 1500;
setup.derivatives.supplier ='sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.method = 'RPM-Differentiation';
setup.displaylevel = 1;
%setup.mesh.phase.fraction = 0.05*ones(1,20);
%%
%---------------------- Solve Problem Using GPOPS2 -----------------------%
tic
output = gpops2(setup);
et = toc;

solution = output.result.solution;
ts = solution.phase.time;
xc = solution.phase.state(:,1:end);
U = solution.phase.control;
%% Plotting -----------------------------
% str1 = sprintf('Koopman simulation with %1.0f sim*%1.0f traj',Nsim,Ntraj);
lw=2;
figure
nexttile
plot(ts,U(:),'linewidth',lw); hold on
ylabel('torque input ', 'interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(ts,xc(:,1),'linewidth',lw)
% title(str1)
ylabel('$u$','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(ts,xc(:,2),'linewidth',lw);
ylabel('$\theta$ yaw','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(ts,xc(:,3),'linewidth',lw);
ylabel('$\dot{\theta}$ yaw rate','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(ts,xc(:,4),'linewidth',lw);hold on
plot(ts,0*xc(:,4),'--r','linewidth',lw);
ylabel('$\psi$ roll','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(ts,xc(:,5),'linewidth',lw);
ylabel('$\dot{\psi}$ roll rate','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

figure
plot(xc(:,1),xc(:,3),'linewidth',lw);
ylabel('$\dot{\theta}$ roll rate','interpreter','latex');
xlabel('u','interpreter','latex');
set(gca,'fontsize',20)


%%
% BEGIN: fContinuous - dynamics(E.O.M), Integrand(power), Path constraints(obstacles) %
%--------------------------------------------------------------------------------------%
function output = fContinuous(input)
global b h M C_s g C_w C_p Ix Iy Iz
const = [Ix;Iy;Iz;h;b;M;g;C_s;C_p;C_w];
f_xu = @eom_grnd_roll_CS; % f(x,u) -- non-linear dynamics with 1 input torque

T1 = input.phase.control(:,1);
u = input.phase.state(:,1); theta = input.phase.state(:,2);
d_theta = input.phase.state(:,3); psi = input.phase.state(:,4);
d_psi = input.phase.state(:,5);
x = [u.';theta.';d_theta.';psi.';d_psi.'];
t = input.phase.time;
ds = f_xu(t.',x,const,T1.');
dq = ds.';
power = (psi.^2) ;%+ 0*u.^2;
% dq=zeros(size(t,1),size(x,1));
% fContinuous Setup
output.integrand = power;
output.dynamics = dq;
% output.path = 0; % c<= 0
end

%%
% BEGIN: fEndpoint - objective(t_final, integral(Energy), Eventgroup(destination) %
%---------------------------------%
function output = fEndpoint(input)

% t_final = input.phase.finaltime;
Energy = input.phase.integral;
% F_pos(1,1) = input.phase.finalstate(:,1); F_pos(2,1) = input.phase.finalstate(:,3);
% F_pos(3,1) = input.phase.finalstate(:,5);
% terminal_cost = norm(F_pos-Xdes);

% output.eventgroup.event = 0;%terminal_cost;
output.objective =  Energy;
end
