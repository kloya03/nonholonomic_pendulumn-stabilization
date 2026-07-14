clc;
clear;

%% Simulation Comparison -----------------
load('koopman_roll_CS.mat') %
% load('koopman1.mat')
tsim = 5; % simulation time
tpred = 0.05; % prediction horizon
rn = 2;
tcont = tpred;%0.1; % control horizon
tc = 0; % current time
nsim = tsim/dt;


A1 = -50; A2 = 32; W1= 20; W2 = 12;
in4 = [A1;W1;A2;W2]; % ode45 input
u_dt = @(i)(A1*sin((W1*i*dt))+A2*cos(W2*i*dt)); % discrete control signal

% u_dt = @(i)(-2*sin((10*i*dt))); % control signal
% inps = [-2;10];
nc_count =[];
% Initial condition
x0 = [0;0;0.1;0.5;0.75];
x_true = x0; Udt = [];
x_truec = x_true;
% Lifted initial condition

xlift0 = double(subs(zm,y,x0));
xliftc = xlift0;
xlift = xlift0;

while tc < tsim
    if tsim-tc <=tcont
        Np = round((tsim-tc)/dt);
        Nc = Np;
    elseif (tcont < tsim-tc)&& (tsim-tc <= tpred)
        Np = round((tsim-tc)/dt);
        Nc = round(tcont/dt);
    else
        Np = round(tpred/dt);
        Nc = round(tcont/dt);
    end

    ix = round(tc/dt):round(tc/dt)+Nc;
    nc_count = [nc_count Nc];
    for j = ix(1,2):ix(1,end)
        Udt = [Udt u_dt(j)];
        xlift = [xlift, Alift*xlift0 + Blift*u_dt(j)]; % Lifted dynamics
        x_true = [x_true, f_ud(0,x_true(:,end),u_dt(j))];
        xlift0 = xlift(:,end);
    end
    xlift0 =  double(subs(zm,y,x_true(:,end)));
    tc = round(tc + Nc*dt,rn)

end
    % open loop Simulate
%     for i = 0:nsim-1
% 
%            % Koopman predictor
%         xlift = [xlift, Alift*xlift(:,end) + Blift*u_dt(i)]; % Lifted dynamics
% 
%            % True dynamics
%         x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
%     end

% closed loop Simulate
for i = 0:nsim-1
    % Koopman predictor
    x_liftt = double(subs(zm,y,x_truec(:,end)));
    xliftc = [xliftc, Alift*x_liftt + Blift*u_dt(i)]; % Lifted dynamics

    % True dynamics
    x_truec = [x_truec, f_ud(0,x_truec(:,end),u_dt(i)) ];
end

x_koop = Clift * xlift; % Koopman predictions
% x_koopc = Clift*xliftc;
[t1,sts1] = ode45(@(t,sts)eom_grnd_roll_CS_T(t,sts,const,in4),[0 tsim],x0);

%% Plotting -----------------------------
% str1 = sprintf('Koopman simulation with %1.0f sim*%1.0f traj',Nsim,Ntraj);
lw=2;
figure
nexttile
plot(t1,(A1*sin((W1*t1))+ A2*cos(W2*t1)),'k','linewidth',lw);  hold on
plot(0:dt:tsim-dt,Udt,'--r','linewidth',lw);
plot(0:dt:tsim-dt,u_dt(0:nsim-1),'--g','linewidth',lw);
ylabel('torque input ', 'interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

figure
nexttile
plot(t1,sts1(:,1),'k','linewidth',lw); hold on
plot([0:nsim]*dt,x_true(1,:),'--b','linewidth',lw); hold on
plot([0:nsim]*dt,x_koop(1,:), '--r','linewidth',lw)
% plot([0:nsim]*dt,x_koopc(1,:), '--g','linewidth',lw)
% title(str1)
ylabel('$u$','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('ODE45','True RK-4','Koopman','location','south');
set(LEG,'interpreter','latex')

nexttile
plot(t1,sts1(:,2),'k','linewidth',lw); hold on
plot([0:nsim]*dt,x_true(2,:),'-b','linewidth',lw); hold on
plot([0:nsim]*dt,x_koop(2,:), '--r','linewidth',lw); hold on
% plot([0:nsim]*dt,x_koopc(3,:), '--g','linewidth',lw)
ylabel('$\theta$ roll ','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(t1,sts1(:,3),'k','linewidth',lw); hold on
plot([0:nsim]*dt,x_true(3,:),'--b','linewidth',lw); hold on
plot([0:nsim]*dt,x_koop(3,:), '--r','linewidth',lw)
% plot([0:nsim]*dt,x_koopc(4,:), '--g','linewidth',lw)
ylabel('$\dot{\theta}$ roll rate','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(t1,sts1(:,4),'k','linewidth',lw); hold on
plot([0:nsim]*dt,x_true(4,:),'-b','linewidth',lw); hold on
plot([0:nsim]*dt,x_koop(4,:), '--r','linewidth',lw); hold on
% plot([0:nsim]*dt,x_koopc(3,:), '--g','linewidth',lw)
ylabel('$\psi$ roll ','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(t1,sts1(:,5),'k','linewidth',lw); hold on
plot([0:nsim]*dt,x_true(5,:),'--b','linewidth',lw); hold on
plot([0:nsim]*dt,x_koop(5,:), '--r','linewidth',lw)
% plot([0:nsim]*dt,x_koopc(5,:), '--g','linewidth',lw)
ylabel('$\dot{\psi}$ roll rate','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
% LEG = legend('True','Koopman','location','southwest');
% set(LEG,'interpreter','latex')
% saveas(figure(1),'Cart_pend_Predictor_comparison')
