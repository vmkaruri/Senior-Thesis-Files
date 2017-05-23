function [x y z]=lorenzo(iter,dt);
x(1)=2.9; y(1)=-1.3; z(1)=25.9;
for i=2:iter
    [x(i) y(i) z(i)]=nextit(x(i-1),y(i-1),z(i-1),dt);
end




function [XO YO ZO]=nextit(X,Y,Z,dt)
    r=28; s=8/3; t=10;
    x1=X+t*(Y-X)*dt/2;
	y1=Y+(X*(r-Z)-Y)*dt/2;
	z1=Z+(X*Y-s*Z)*dt/2;
	XO=X+t*(y1-x1)*dt;
	YO=Y+(x1*(r-z1)-y1)*dt;
	ZO=Z+(x1*y1-s*z1)*dt;