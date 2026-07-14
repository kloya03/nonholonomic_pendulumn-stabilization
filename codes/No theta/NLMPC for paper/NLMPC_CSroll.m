clear all;
clc;
tic
nx = 4;
ny = 4;
nu = 1;
% load('koopman_roll_CS_not_randn12.mat')
f_xu1 = @eom_grnd_roll_CS_red1; % f(x,u) -- non-linear dynamics with 1 input torque
global b h M C_s g C_w C_p I_xx I_zz I_yy
b = 0.15; %m
h = 0.05; % m
M = 2.5; %kg
C_s = 0.85;
g =9.8;
C_w = 0.85;
C_p = 0.9;

% Moment of Inertia ellipsoid
ae = b; be = b/1.5; ce = be;
I_zz = 0.2*M*(ae^2 + be^2);
I_xx = 0.2*M*(ce^2 + be^2) + M*h^2;
I_yy = 0.2*M*(ce^2 + ae^2) + M*h^2;

%% nlobj
nlobj = nlmpc(nx,ny,nu);
Ts = 0.01;
np = 50;
nc =20;
nlobj.Ts = Ts;
nlobj.PredictionHorizon = np;
nlobj.ControlHorizon = nc;

%% Model
nlobj.Model.StateFcn = @(x,u)eom_grnd_roll_CS_red1(x,u);
nlobj.Model.IsContinuousTime = true;
% nlobj.Model.NumberOfParameters = 1;

%% Constraints
for ct = 1:nx
    nlobj.states(ct).Min = -100; %% also can include ratemin max
    nlobj.states(ct).Max = 100;
end
nlobj.MV.Min = -75;
nlobj.MV.Max = 75;

%% Cost function
nlobj.Optimization.CustomCostFcn = @Nl_costfcn;
nlobj.Optimization.ReplaceStandardCost = true;

%% initial state and inputs
x0 = [0;0;-0.1;0];
u0 = zeros(nu,1);
validateFcns(nlobj,x0,u0);

%% Solve for N-L MPC
nloptions = nlmpcmoveopt;
nloptions.Parameters = {Ts};
optns = optimset('display','iter','MaxIter', 500, 'MaxFunEvals', 100000);
optim.options.Fmincon = optns;
Duration = 2;
xHistory = x0;
mv = u0;
mu = u0;
% x = x0;
% for ct = 1:(Duration/Ts)
%     ct
%     % Compute optimal control moves
%     [mv,nloptions,info] = nlmpcmove(nlobj,x,mv,[],[],nloptions);
%     % Implement first optimal control move
%     x = x + Ts*f_xu(x,mv);
%     % save control actions
%     mu = [mu mv];
%     % Save plant states
%     xHistory = [xHistory x];
% end
xc = x0;
for ct = 1:(Duration/(nc*Ts))
    ct
    % Compute optimal control moves
    [mv,nloptions,info] = nlmpcmove(nlobj,xc(:,end),mv,[],[],nloptions);
    % Implement first optimal control move
    uc = nloptions.MV0(1:nc,1);
    for i = 1:nc
    xc = [xc xc(:,end) + Ts*f_xu1(xc(:,end),uc(i,1))];
    end
    % save control actions
    mu = [mu;uc];
    % Save plant states
%     xHistory = [xHistory x];
end
 xHistory = xc;
% [~,~,info] = nlmpcmove(nlobj,x0,u0);
et = toc;

figure
plot(0:Ts:Duration,mu)
figure
subplot(2,2,1)
plot(0:Ts:Duration,xHistory(1,:))
xlabel('time')
ylabel('u')
title('velocity')
subplot(2,2,2)
plot(0:Ts:Duration,xHistory(2,:))
xlabel('time')
ylabel('thetadot')
title('angular velovity')
subplot(2,2,3)
plot(0:Ts:Duration,xHistory(3,:))
xlabel('time')
ylabel('psi')
title('roll angle')
subplot(2,2,4)
plot(0:Ts:Duration,xHistory(4,:))
xlabel('time')
ylabel('psidot')
title('roll angle velocity')

%%
% dt = Ts;
% tsim = duration;
% xRk = x0;
% for i = 1:tsim/dt
%     xRk = [xRk, f_ud(0,xRk(:,end),mu(i))];
% end
% figure
% nexttile
% plot(0:dt:tsim-dt,mu(:));hold on;
% plot(0:dt:tsim-dt,T3(:))
% ylabel('total control')
% nexttile
% plot(0:dt:tsim,xRk(1,:))
% ylabel('u')
% title('RK4 check')
% nexttile
% plot(0:dt:tsim,xRk(2,:))
% ylabel('$\dot{\theta}$','interpreter','latex')
% % nexttile
% % plot(xRk(1,:),xRk(3,:))
% nexttile
% plot(0:dt:tsim,xRk(3,:))
% ylabel('$\psi$','interpreter','latex')
% nexttile
% plot(0:dt:tsim,xRk(4,:))
% ylabel('$\dot{\psi}$','interpreter','latex')


%% Ode45 plots
const(:,1) =  [I_xx;I_yy;I_zz;h;b;M;g;C_s;C_p;C_w];
[t1,sts1] = ode45(@(t,sts)eom_grnd_roll_CS_U(t,sts,const,mu),[0 Duration],x0);
Uode = [];
for i = 1:size(t1,1)
    k = t1(i);
    k = ceil(1000*t1(i));
    if k==0
        Uo = mu(1);
    else
        Uo = mu(k);
    end
    Uode = [Uode;Uo];
end
figure
nexttile
plot(t1,Uode(:));hold on
% plot(t1,A*sin(Omega*t1),'r')
ylabel('total control')
nexttile
plot(t1,sts1(:,1))
ylabel('u')
title('ode45 check')
nexttile
plot(t1,sts1(:,2))
ylabel('$\dot{\theta}$','interpreter','latex')
% nexttile
% plot(sts2(:,1),sts2(:,3))
nexttile
plot(t1,sts1(:,3))
ylabel('$\psi$','interpreter','latex')
nexttile
plot(t1,sts1(:,4))
ylabel('$\dot{\psi}$','interpreter','latex')

%%
function J = Nl_costfcn(X,U,e,data)
vel = norm(X(:,1)-2);
stability = norm(X(:,3));
costt = norm(U(:,data.MVIndex(1)));
J = stability + vel + 1e-06*costt;
end