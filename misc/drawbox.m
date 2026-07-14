function drawbox(Xo,Yo,Zo,theta, phi)

l = 2; 
w = 1;
h = 1;

[X,Y,Z] = ellipsoid(0,0,0,0.15,0.07,0.07);
tt = linspace(0, 2*pi);
r = 0.07;
x = -0.2+r*cos(tt);
z = 0+r*sin(tt);
y=0;
% surf(X,Y,Z)
% X = [l,l,l,l;           % front
%     -l, -l, -l, -l;     % back
%     l, -l, -l, l;       % top
%     l, -l, -l, l;       % bottom
%     l , l, -l, -l;      % right
%     l , l, -l, -l];     % left
% Y = [w,w,-w,-w;
%      w,w,-w,-w;
%      w, w, -w, -w;
%      w, w, -w, -w;
%      w, w, w, w;
%      -w, -w, -w, -w;];
% Z = [h,-h,-h,h;
%      h,-h,-h,h; 
%      h, h, h, h;
%      -h, -h, -h, -h;
%      h, -h, -h, h;
%      h, -h, -h, h];


Rz = [cos(theta) , -sin(theta) , 0; 
      sin(theta) ,cos(theta), 0; 
      0 , 0 , 1];
  
Rx = [1 , 0 , 0; 
      0 , cos(phi), -sin(phi) ;
      0 , sin(phi), cos(phi)];
  
R = Rx*Rz;

Xs = Xo + R(1,1)*X + R(1,2)*Y + R(1,3)*Z;
Ys = Yo + R(2,1)*X + R(2,2)*Y + R(2,3)*Z;
Zs = Zo + R(3,1)*X + R(3,2)*Y + R(3,3)*Z;

xs = Xo + R(1,1)*x + R(1,2)*y + R(1,3)*z;
ys = Yo + R(2,1)*x + R(2,2)*y + R(2,3)*z;
zs = Zo + R(3,1)*x + R(3,2)*y + R(3,3)*z;

surf(Xs(1:13,:),Ys(1:13,:),Zs(1:13,:));
patch(Xs(11,:),Ys(11,:),Zs(11,:),0.23+0*Zs(11,:));
patch(xs, ys, zs, 'k')
% fill3(Xs',Ys',Zs','blue','FaceAlpha',0.5);
% surf([Xs;Xs'],[Ys;Ys'],[Zs;Zs'],'blue','FaceAlpha',0.5);
% axis([-5 5 -5 5 -5 5])
end
% [X,Y,Z] = ellipsoid(0,0,0,1.5,0.75,0.75);
% surf(X(9:13,:),Y(9:13,:),Z(9:13,:))
% hold on
% patch(X(13,:),Y(13,:),Z(13,:),Z(13,:));
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
% for i=1:21
%     plot3(X(13,i),Y(13,i),Z(13,i),'*r')
%     hold on
% end