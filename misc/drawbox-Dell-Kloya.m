function drawbox(Xo,Yo,Zo,theta, phi)

l = 2; 
w = 1;
h = 0.2;

X = [l,l,l,l;           % front
    -l, -l, -l, -l;     % back
    l, -l, -l, l;       % top
    l, -l, -l, l;       % bottom
    l , l, -l, -l;      % right
    l , l, -l, -l];     % left
Y = [w,w,-w,-w;
     w,w,-w,-w;
     w, w, -w, -w;
     w, w, -w, -w;
     w, w, w, w;
     -w, -w, -w, -w;];
Z = [h,-h,-h,h;
     h,-h,-h,h; 
     h, h, h, h;
     -h, -h, -h, -h;
     h, -h, -h, h;
     h, -h, -h, h];


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

fill3(Xs',Ys',Zs','blue','FaceAlpha',0.5);

end

