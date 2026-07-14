clc;
clear;

%% Simulation Comparison -----------------
load('koopman_roll_CS_cossine.mat') %
% load('koopman1.mat')
tsim = 3; % simulation time
tpred = 0.01; % prediction horizon
rn = 2;
tcont = tpred;%0.1; % control horizon
tc = 0; % current time
nsim = tsim/dt;

u_dt = @(i)(-2*sin((10*i*dt))); % control signal
inps = [-2;10];
nc_count =[];
% Initial condition
x0 = [0;0;0.1;0.15;0];
x0e = [0;0;0.1;0.15;0;cos(0.15);sin(0.15)];
x_true = x0; Udt = [];
x_truec = x_true;
% Lifted initial condition

xlift0 = double(subs(zm,y,x0e));
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
    xe_true = [x_true;cos(x_true(4,:));sin(x_true(4,:))];
    xlift0 =  double(subs(zm,y,xe_true(:,end)));
    tc = round(tc + Nc*dt,rn)

end
%     % open loop Simulate
%     for i = 0:nsim-1
%
%         %     Koopman predictor
%         xlift = [xlift, Alift*xlift(:,end) + Blift*u_dt(i)]; % Lifted dynamics
%
%         %     True dynamics
%         x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
%     end

% % closed loop Simulate
% for i = 0:nsim-1
%     % Koopman predictor
%     x_liftt = double(subs(zm,y,x_truec(:,end)));
%     xliftc = [xliftc, Alift*x_liftt + Blift*u_dt(i)]; % Lifted dynamics
% 
%     % True dynamics
%     x_truec = [x_truec, f_ud(0,x_truec(:,end),u_dt(i)) ];
% end

x_koop = Clift * xlift; % Koopman predictions
% x_koopc = Clift*xliftc;
[t1,sts1] = ode45(@(t,sts)eom_grnd_roll_CS_ode45(t,sts,const,inps),[0 tsim],x0);

%% Plotting -----------------------------
% str1 = sprintf('Koopman simulation with %1.0f sim*%1.0f traj',Nsim,Ntraj);
lw=2;
figure
nexttile
plot(t1,inps(1,1)*sin(inps(2,1)*t1),'k','linewidth',lw);  hold on
plot(0:dt:tsim-dt,Udt,'--r','linewidth',lw);
% plot(0:dt:tsim-dt,u_dt(0:nsim-1),'--g','linewidth',lw);
ylabel('torque input ', 'interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

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
% plot([0:nsim]*dt,x_koopc(4,:), '--g','linewidth',lw)
ylabel('$\dot{\psi}$ roll rate','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
% LEG = legend('True','Koopman','location','southwest');
% set(LEG,'interpreter','latex')
% saveas(figure(1),'Cart_pend_Predictor_comparison')
