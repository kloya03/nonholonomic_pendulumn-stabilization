%% %%%%%%%%% MODEL PREDICTIVE CONTROL %%%%%%%%% %%
clc;
clear;
%% Initialize
load('koopman_roll_CS_not_norm.mat','Nlift','Alift','Blift','dt','n','m','x','zm','f_ud','f_xu','const')
% Alift, Blift, Clift, dt, umin, umax, yrr, n, m, Nlift
% xs =  state
% x(1) = lifted_fun(xs)
% X = [x(1) x(2) x(3) ...x(Np) u(1) u(2) u(3)...u(Np-1)]

tsim = 2; % simulation time
tpred = 0.01; % prediction horizon
tcont = 0.005;%0.05; % control horizon
tc = 0; % current time
ur = 7.758;
x0 = [0;0;0.35;0.15]; % initial state
xlift0 = double(subs(zm,x,x0)); % initial lifted state
xc = x0; U = [];nsim = 0;
eps = 1e-10;
fal = [];     T3 = []; Ut=[];
A= 0; Omega = 15; % for uc = 4;
maxA=50;
um=maxA-A;

%% Limit cycle plot
% [t2,sts2] = ode45(@(t,sts)eom_CS(t,sts,const,A*sin(Omega*t)),[0 30],[2;0;0.2;0.35;0.1]);
% figure
% nexttile
% plot(t2,A*sin(Omega*t2))
% ylabel('T1')
% title('A*sin(\Omega*t) effect')
% nexttile
% plot(t2,sts2(:,3))
% ylabel('u')
% nexttile
% plot(t2,sts2(:,5))
% ylabel('$\dot{\theta}$','interpreter','latex')
% nexttile
% plot(sts2(:,3),sts2(:,5))

%% MPC Begins
tic
while tc < tsim -eps
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
    wd =0*rand(0,1);
    T2 = [];
    %% weights on states
    %qs = [0;0;0;1;0]; ql = double(subs(zm,y,qs)); %
    eps1 = 10/(ur^2);
    ql = zeros(Nlift,1);
    ql(4,1) =1;
    ql(1,1) = 100;
    %     ql(1,1) = eps1*2*ur^2;
    %     ql(2,1) = eps1*2; %q(7,1) = 1;
    R =1e-07 ;%0*0.0001*(0.1)/um^2;
    Qx = sparse(repmat(ql,Np,1)); Qu = sparse(repmat(R,Np,1));
    QNp = (diag([Qx;Qu]));
    
    f = zeros(Np*(Nlift+m),1);
    %     for ff = 1:Np
    %         %         f((ff-1)*Nlift+4,1) = 10;
    % %         f((ff-1)*Nlift+2,1) = -eps1*2*ur;
    %         %         f((ff-1)*Nlift+6,1) = eps1;
    %     end
    %% Constraint
    beq1=[];
    Aeq1 = sparse(zeros(Nlift*(Np),(Nlift+m)*Np));
    Aeqc = mat2cell(Aeq1,Nlift*ones(1,Np),[Nlift*ones(1,Np) m*ones(1,Np)]);
    for i=1:Np
        if i<Np
            Aeqc{i+1,i} = Alift;
        end
        Aeqc{i,i} = -eye(Nlift,Nlift);
        Aeqc{i,Np+i} = Blift;
        Ti = A*sin(Omega*(tc+(i-1)*dt));
        T2 = [T2; Ti];
        beq1 = [beq1;Blift*Ti];
    end
    Aeq = cell2mat(Aeqc);
    Beq = [-Alift;zeros((Nlift)*(Np-1),Nlift)]*xlift0 - beq1;
    
    %% Quadprog setup
    lb = [-inf*ones(Nlift*Np,1);-um*ones(Np,1)];
    ub = [inf*ones(Nlift*Np,1);um*ones(Np,1)];
    %     options = optimoptions('quadprog','Algorithm','active-set');
    %     ini_X = xlift0;
    [X,fval,exitflag,output] = quadprog(QNp,f,[],[],Aeq,Beq,lb,ub);%,ini_X,options);
    
    %% Trajectory Extraction
    fal = [fal fval];
    %     psii = [];
    %     for kk = 1:tsim/dt
    %     psii = [psii;X(Nlift*(kk-1)+4,1)];
    %     end
    
    u_koop = X(Np*Nlift +1:Np*Nlift+m*Nc,1);
    
    for ix = 1:Nc
        T3 = [T3 T2(ix)];
        U = [U u_koop(ix)];
        xc = [xc, f_ud(0,xc(:,end),T2(ix)+u_koop(ix))];
        %         xc = [xc, f_ud(0,xc(:,end),u_koop(ix))];
    end
    xlift0 =  double(subs(zm,x,xc(:,end)));
    nsim = nsim+Nc;
    tc = round(tc + Nc*dt,3)
end
et=toc;
% tt = dt:dt:tsim;
% cntrl_MPC = [tt.' U(:)];

%% Runge -Kutta plots
xRk = x0;
for i = 1:tsim/dt
    xRk = [xRk, f_ud(0,xRk(:,end),U(i)+T3(i))];
end
figure
nexttile
plot(0:dt:tsim-dt,U(:));hold on;
plot(0:dt:tsim-dt,T3(:))
ylabel('total control')
nexttile
plot(0:dt:tsim,xRk(1,:))
ylabel('u')
title('RK4 check')
nexttile
plot(0:dt:tsim,xRk(2,:))
ylabel('$\dot{\theta}$','interpreter','latex')
% nexttile
% plot(xRk(1,:),xRk(3,:))
nexttile
plot(0:dt:tsim,xRk(3,:))
ylabel('$\psi$','interpreter','latex')
nexttile
plot(0:dt:tsim,xRk(4,:))
ylabel('$\dot{\psi}$','interpreter','latex')

%% Ode45 plots
[t1,sts1] = ode45(@(t,sts)eom_grnd_roll_CS_U(t,sts,const,U(:),A,Omega),[0 tsim],x0);
Uode = [];
for i = 1:size(t1,1)
    k = t1(i);
    k = ceil(1000*t1(i));
    if k==0
        Uo = U(1);
    else
        Uo = U(k);
    end
    Uode = [Uode;Uo];
end
figure
nexttile
plot(t1,Uode(:));hold on
plot(t1,A*sin(Omega*t1),'r')
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


%% Plotting -----------------------------
% str1 = sprintf('Koopman simulation with %1.0f sim*%1.0f traj',Nsim,Ntraj);
tsimm = tsim-0.05;
lw=4;
lw1 = 35;
lw2=60;
figure%('Units','normalized','Position',[0 0 1 1])
tiledlayout(1,2)
nexttile
plot(0:dt:tsim-dt,U(:),'linewidth',lw); hold on
% plot(0:dt:tsim-dt,T3,'r','linewidth',lw);hold on
% plot(0:dt:tsim-dt,T2,'r','linewidth',lw/2); hold on
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',lw1)
ylabel('$\tau$', 'interpreter','latex','fontsize',lw2);
axis([0-0.02 tsimm -inf inf])

% figure('Units','normalized','Position',[0 0 1 1])
% tiledlayout(1,2)
nexttile %% psi
% plot(t1,sts1(:,4),'k','linewidth',lw); hold on
plot(0:dt:tsim,xc(3,:),'linewidth',lw);hold on
plot(0:dt:tsim,0*xc(3,:),'--r','linewidth',lw);
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',lw1)
ylabel('$\psi$','interpreter','latex','fontsize',lw2);
axis([0-0.02 tsimm -0.01 inf])

figure%('Units','normalized','Position',[0 0 1 1])
tiledlayout(1,3)
nexttile %% u
% plot(t1,sts1(:,1),'k','linewidth',lw); hold on
plot(0:dt:tsim,xc(1,:),'linewidth',lw)
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',lw1)
ylabel('$u$','interpreter','latex','fontsize',lw2);
axis([0-0.04 tsimm -inf inf])

% nexttile
% % plot(t1,sts1(:,2),'k','linewidth',lw); hold on
% plot(0:dt:tsim,xc(2,:),'linewidth',lw);
% ylabel('$\theta$ yaw','interpreter','latex');
% xlabel('Time [s]','interpreter','latex');
% set(gca,'fontsize',20)

% figure('Units','normalized','Position',[0 0 1 1])
% tiledlayout(1,2)
nexttile  %% theta_dot
% plot(t1,sts1(:,3),'k','linewidth',lw); hold on
plot(0:dt:tsim,xc(2,:),'linewidth',lw);
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',lw1)
ylabel('$\dot{\theta}$','interpreter','latex','fontsize',lw2);
axis([0-0.04 tsimm -1 inf])

nexttile  %% psi_dot
% plot(t1,sts1(:,5),'k','linewidth',lw); hold on
plot(0:dt:tsim,xc(4,:),'linewidth',lw);
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',lw1)
ylabel('$\dot{\psi}$','interpreter','latex','fontsize',lw2);
axis([0-0.04 tsimm -inf 1])


% figure
% nexttile
% plot(xc(1,:),xc(3,:),'linewidth',lw);
% ylabel('$\dot{\theta}$ roll rate','interpreter','latex');
% xlabel('u','interpreter','latex');
% set(gca,'fontsize',20)