function [a1,a2]=perturbation(x,y,z,Re,mu)
%地球引力加速度（含非球形摄动加速度） 
r=sqrt(x^2+y^2+z^2);
[e1,e2,e3]=transformation(x,y,z);
p=asin(z/r);
if(y>0)
    q=acos(x/sqrt(x^2+y^2));
else
    q=2*pi-acos(x/sqrt(x^2+y^2));
end
a1=-mu*e1/r^2;
a2=0;
for n=2:4
    b=(Re/r)^n;
    for m=0:n
        [C,S]=CS(n,m);
        c1=(1+n)*cos(p)*P(n,m,sin(p))*(C*cos(m*q)+S*sin(m*q))*e1;
        c2=m*P(n,m,sin(p))*(C*sin(m*q)-S*cos(m*q))*e2;
        c3=(n*sin(p)*P(n,m,sin(p))-sqrt((2*n+1)*(n+m)*(n-m)/(2*n-1))*P(n-1,m,sin(p)))*(C*cos(m*q)+S*sin(m*q))*e3;
        c=c1+c2+c3;
        a2=a2-mu*b*c/(r^2*cos(p));
    end
end
end