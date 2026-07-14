%% %%%%%%%%% MODEL PREDICTIVE CONTROL %%%%%%%%% %%
clc;
clear;
load('koopman1.mat')
% Alift, Blift, Clift, dt, umin, umax, yrr, n, m, Nlift
% xs =  state
% x(1) = lifted_fun(xs)
% X = [x(1) x(2) x(3) ...x(Np) u(1) u(2) u(3)...u(Np-1)]

tsim = 2; % simulation time
tpred = 0.1; % prediction horizon
tcont = 0.01; % control horizon
tc = 0; % current time
% r = 1;
% nsim = (tsim)/dt;
% Np = round(tpred/dt);
% Nc = round(tcont/dt);

x0 = [0;0;1.9;1.12;0.9]; % initial state
xlift0 = double(subs(zm,y,x0)); % initial lifted state
xc = x0; U = [];nsim = 0;
eps = 1e-10;
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
    % weights on states
    %         qs = [0;0;0;1;0]; ql = double(subs(zm,y,qs)); %ql(1,1) = 0;
    ql = zeros(Nlift,1); ql(5,1) = 1; %q(2,1) = 1;%r^2; q(2,1) = -2*r; q(7,1) = 1; 
    R = 0; Qx = sparse(repmat(ql,Np,1)); Qu = repmat(R,Np,1);
    QNp = sparse(diag([Qx;Qu]));
    
    % Constraint
    Aeq1 = sparse(zeros(Nlift*(Np),(Nlift+m)*Np));
    Aeqc = mat2cell(Aeq1,Nlift*ones(1,Np),[Nlift*ones(1,Np) m*ones(1,Np)]);
%     syms a b I
%     Np = 5;
%     Alift = a;
%     Blift = b;
%     Nlift = 1;
    for i=1:Np
        if i<Np
            Aeqc{i+1,i} = Alift;
        end
        Aeqc{i,i} = -eye(Nlift,Nlift);
        Aeqc{i,Np+i} = Blift;
    end
    Aeq = cell2mat(Aeqc);
    Beq = [-Alift;zeros((Nlift)*(Np-1),Nlift)]*xlift0;
    
    lb = [-1000000*ones(Nlift*Np,1);-75*ones(Np,1)];
    ub = [1000000*ones(Nlift*Np,1);75*ones(Np,1)];
%     options = optimoptions('quadprog','Algorithm','active-set');
%     ini_X = xlift0;
    [X,fval,exitflag,output] = quadprog(QNp,[],[],[],Aeq,Beq,lb,ub);%,ini_X,options);
    
    u_koop = X(Np*Nlift +1:Np*Nlift+m*Nc,1);
    for ix = 1:Nc
        xc = [xc, f_ud(0,xc(:,end),u_koop(ix))];
    end
    xlift0 =  double(subs(zm,y,xc(:,end)));
    U = [U u_koop.'];
    nsim = nsim+Nc;
    tc = round(tc + Nc*dt,2)
end
et=toc
%% Plotting -----------------------------
% str1 = sprintf('Koopman simulation with %1.0f sim*%1.0f traj',Nsim,Ntraj);
lw=2;
figure
nexttile
plot(0:dt:tsim,[zeros(m,1);U(:)],'linewidth',lw); hold on
ylabel('torque input ', 'interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(0:dt:tsim,xc(1,:),'linewidth',lw)
% title(str1)
ylabel('$u$','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(0:dt:tsim,xc(2,:),'linewidth',lw);
ylabel('$\theta$ yaw','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(0:dt:tsim,xc(3,:),'linewidth',lw);
ylabel('$\dot{\theta}$ yaw rate','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(0:dt:tsim,xc(4,:),'linewidth',lw);hold on
plot(0:dt:tsim,0*xc(4,:),'--r','linewidth',lw);
ylabel('$\psi$ roll','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

nexttile
plot(0:dt:tsim,xc(5,:),'linewidth',lw);
ylabel('$\dot{\psi}$ roll rate','interpreter','latex');
xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

figure
plot(xc(1,:),xc(3,:),'linewidth',lw);
ylabel('$\dot{\theta}$ roll rate','interpreter','latex');
xlabel('u','interpreter','latex');
set(gca,'fontsize',20)


% saveas(figure(1),'Predictor_comparison')