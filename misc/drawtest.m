clear
clc
close all
load("Paper_result_disturb_005.mat",'xc')
%%
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
const(:,1) = [Ix;Iy;Iz;h;b;M;g;C_s;C_p;C_w];
% sts_c =[s(1,1);s(2,1);q_dot(1,1);s(3,1);q_dot(3,1)];
% A = 2; Omega = 0.5; U = 0.8;
% inp_W = [2;15;0;0];
x0 = [0;0;0;0;];%0;0;0];
tspan = 0:0.001:2;
[t1,sts2] = ode45(@(t,sts)eom_grnd_roll_CS_U(t,sts,const,1),tspan,x0);
sts1(:,3) = sts2(:,1); sts1(:,5) = sts2(:,2); sts1(:,6) = sts2(:,3);

dt = 0.001;
% sts1(:,3) = xc(1,:).'; sts1(:,5) = xc(2,:).';
% sts1(:,6) = xc(3,:).';
sts1(1,1) = x0(1,1); sts1(1,2) = x0(2,1);
sts1(1,4) = x0(4,1);
for i=1:size(tspan,2)-1
    sts1(i+1,4) = sts1(i,4) + sts2(i+1,2)*dt;
    sts1(i+1,1) = sts1(i,1) + (sts1(i,1)*cos(sts1(i,4)) - b*sts1(i,4)*sin(sts1(i,4)))*dt;
    sts1(i+1,2) = sts1(i,2) + (sts1(i,1)*sin(sts1(i,4)) + b*sts1(i,4)*cos(sts1(i,4)))*dt;
end
% for i=1:size(tspan,2)-1
%     sts1(i+1,4) = sts1(i,4) + xc(2,i+1)*dt;
%     sts1(i+1,1) = sts1(i,1) + (xc(1,i)*cos(sts1(i,4)) - b*sts1(i,4)*sin(sts1(i,4)))*dt;
%     sts1(i+1,2) = sts1(i,2) + (xc(1,i)*sin(sts1(i,4)) + b*sts1(i,4)*cos(sts1(i,4)))*dt;
% end
%%

f1 = figure('color','w');

% nt = 100;
% t = linspace(0,10,nt);
t = tspan;
nt = size(t,2);

X = sts1(:,1);
Y = sts1(:,2);
Z = zeros(size(X));
Th = sts1(:,4);
Phi = sts1(:,6);
fram = nt/tspan(1,end);
% v = VideoWriter('BoxMovie005.mp4','MPEG-4');
% v.Quality = 100;
% v.FrameRate = 160;
% open(v)

for ii = 1:nt
    
    clf
    hold on 
    plot3(X(1:ii), Y(1:ii), Z(1:ii),'-r')
    plot3(X(ii), Y(ii), Z(ii), 'ok','markerfacecolor','k','markersize',2)
    drawbox(X(ii),Y(ii),Z(ii),Th(ii), Phi(ii))
    axis equal 
%     axis([-0.15 1.5 -0.2 0.2 -0.2 0.2])
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on 
    drawnow 
    
%     frame = getframe(f1);
%     im = frame2im(frame);
%     writeVideo(v,frame);
    
end

% close(v)