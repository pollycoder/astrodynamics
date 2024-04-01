function [v1,v2,a,e,theta,iter]=LambSol(r1,r2,t,mu,way,N,branch,MaxIter,epsilon)
%LambSol Lambert问题求解器
%输入参数r1,r2,t,lw,N,branch,MaxIter,mu,epsilon:
%   r1(3), r2(3):初始点和末端点的位置列矢量分量X、Y、Z,单位：米
%   t:转移时间,单位：秒
%   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
%      默认为地心引力系数3.98600441800e+14
%   way:=0表示选择逆时针运行的轨道,否则选择顺时针运行的轨道，默认为逆时针
%   branch:用于多圈转移(N>=1),=1表示选择左支,否则选择右支.默认为左支
%   MaxIter:最大迭代次数,默认取60
%   epsilon:最小精度,默认为1.0e-11
%输出参数[v1,v2,a,e,theta,iter]:
% 	v1(3),v2(3):初始点和末端点的速度列矢量分量,单位：米每秒
% 	a:转移轨道半长轴,单位：米.
%   e:转移轨道偏心率
%   theta:转移角,弧度
%   iter:实际迭代次数
%注:输出参数为一个时,仅输出转移轨道半长轴

% *             清华大学航天航空学院动力学与控制研究室博士生                 *
% *                蒋方华(jiangfh@tsinghua.edu.cn)                      *
% *                 最近修改: 2017.4.4                                    *
% 注：此程序仅对Dario Izzo的程序作了一些修改。对于时间很短的多圈问题，有可能不存
% 在那样的椭圆轨道，而Izzo的程序中并没有给出判断，计算的结果有时是不对的。原程序
% 中对双曲线轨道的半长轴取负值，在此改成正值，并将输出参数p改为偏心率e.
%Programmed by:         Dario Izzo (Advanced Concept Team) Date:
%28/05/2004 Revision:              1 Tested by:             ----------
%please report bugs to dario.izzo@esa.int

%Preliminary control on the function call
if nargin < 9
  epsilon=1.0e-13;
  if nargin<8
      MaxIter=60;
      if nargin<7
          branch=1;
          if nargin<6
              N=0;
              if nargin<5
                  way=0;
                  if nargin<4
                      mu=3.98600441800e+14;
                      if nargin < 3
                          error('ORBIT:LambSol:NotEnoughInputs',...
                              'Not enough input arguments.  See LambSol.'); 
                      end
                  end
              end
          end
      end
  end
end
if t<=0
    warning('ORBIT:LambSol:NotEnoughInputs',...
        'Not NotSuitableInputs inputs.  See LambSol.');
    v1=NaN;
    v2=NaN;
    return
end
if((length(r1)~=length(r2))||(length(r2)>3)||(length(r2)<2))
    error('ORBIT:LambSol:NotSuitableInputs',...
          '两个位置矢量维数应相等,2维或3维.  See LambSol.');
end   
if(MaxIter<1)
    error('ORBIT:LambSol:NotSuitableInputs',...
          '最大迭代次数不能小于1.  See LambSol.');
end

%epsilon=1e-11;  %Increasing the tolerance does not bring any advantage as the 
%precision is usually greater anyway (due to the rectification of the tof
%graph) except near particular cases such as parabolas in which cases a
%lower precision allow for usual convergence.


%Non dimensional units
R=norm(r1);
V=sqrt(mu/R);
T=R/V;

%working with non-dimensional radii and time-of-flight
r1=r1/R;
r2=r2/R;
t=t/T;                     

%Evaluation of the relevant geometry parameters in non dimensional units
r2mod=norm(r2);
theta=real(acos(dot(r1,r2)/r2mod)); %the real command is useful when theta is very 
 %close to pi and the acos function could return complex numbers

crossr1r2=cross(r1,r2);
if(crossr1r2(3)<0.0)
    theta=2.0*pi-theta;
end
if(way)
    theta=2.0*pi-theta;
end
lw=0;
if(theta>pi)
    lw=1;
end

c=sqrt(1+r2mod*r2mod-2.0*r2mod*cos(theta)); %non dimensional chord
s=(1.0+r2mod+c)/2.0;                      %non dimensional semi-perimeter
%am=s/2;                               %minimum energy ellipse semi major axis
lambda=sqrt(r2mod)*cos(theta/2.0)/s;    %lambda parameter defined in BATTIN's book



%We start finding the log(x+1) value of the solution conic:
%%NO MULTI REV --> (1 SOL)
if N==0
    inn1=-0.5233;    %first guess point
    inn2=0.5233;     %second guess point
    x1=log(1.0+inn1);
    x2=log(1.0+inn2);
    y1=log(x2tof(inn1,s,c,lw,N))-log(t);
    y2=log(x2tof(inn2,s,c,lw,N))-log(t);
    
    %Newton iterations
    err=1.0;
    i=0;
    while ((err>epsilon)&& (i<MaxIter) && (y1~=y2))
        i=i+1;
        xnew=(x1*y2-y1*x2)/(y2-y1);
        ynew=log(x2tof(exp(xnew)-1,s,c,lw,N))-log(t);
        x1=x2;
        y1=y2;
        x2=xnew;
        y2=ynew;
        err=abs(x1-xnew);
    end
    if(i>=MaxIter)
        error('ORBIT:LambSol:NoSolutions',...
          '已经达到最大迭代次数，无解.  See LambSol.');       
    else
        iter=i;
        x=exp(xnew)-1.0;
    end    
    
    %%MULTI REV --> (2 SOL) SEPARATING RIGHT AND LEFT BRANCH
else 
    if (branch==1)
        inn1=-0.5234;
        inn2=-0.2234;
    else
        inn1=0.7234;
        inn2=0.5234;
    end
    x1=tan(inn1*pi/2.0);
    x2=tan(inn2*pi/2.0);
    y1=x2tof(inn1,s,c,lw,N)-t;
    
    y2=x2tof(inn2,s,c,lw,N)-t;
    err=1.0;
    i=0;
    
    %Newton Iteration
    while ((err>epsilon) && (i<MaxIter) && (y1~=y2))
        i=i+1;
        xnew=(x1*y2-y1*x2)/(y2-y1);
        ynew=x2tof(atan(xnew)*2.0/pi,s,c,lw,N)-t;
        x1=x2;
        y1=y2;
        x2=xnew;
        y2=ynew;
        err=abs(x1-xnew);	    
    end
    if(i>=MaxIter)
        error('ORBIT:LambSol:NoSolutions',...
          '已经达到最大迭代次数，无解.  See LambSol.');       
    else
        x=atan(xnew)*2.0/pi;
        iter=i;
    end
end
[v1,v2,a,e]=x2sol(x,s,c,lambda,lw,r2mod,theta,r1,r2);
v1=v1*V;
v2=v2*V;
a=a*R;
siz=size(v1);
if(siz(1)<=1)
    v1=v1';
    v2=v2';
end
if nargout==1
    v1=a;
end
%The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
%now need the conic. As for transfer angles near to pi the lagrange
%coefficient technique goes singular (dg approaches a zero/zero that is
%numerically bad) we here use a different technique for those cases. When
%the transfer angle is exactly equal to pi, then the ih unit vector is not
%determined. The remaining equations, though, are still valid.

function [v1,v2,a,e]=x2sol(x,s,c,lambda,lw,r2mod,theta,r1,r2)
am=s/2.0;
a=am/(1-x*x);                       %solution semimajor axis
%calcolo psi
if x<1.0 %ellisse
    beta=2.0*asin(sqrt((s-c)/2.0/a));
    if lw
        beta=-beta;
    end
    alfa=2.0*acos(x);
    psi=(alfa-beta)/2.0;
    eta2=2.0*a*sin(psi)^2/s;
    eta=sqrt(eta2);
else %iperbole
    beta=2.0*asinh(sqrt((c-s)/2.0/a));
    if lw
        beta=-beta;
    end
    alfa=2.0*acosh(x);
    psi=(alfa-beta)/2.0;
    eta2=-2.0*a*sinh(psi)^2/s;
    eta=sqrt(eta2);
end
p=r2mod/am/eta2*sin(theta/2.0)^2;     %parameter of the solution
e=sqrt(1.0-p/a);
sigma1=1.0/eta/sqrt(am)*(2.0*lambda*am-(lambda+x*eta));

ih=cross(r1,r2);
ih=ih./norm(ih);
if lw
    ih=-ih;
end

vr1 = sigma1;
vt1 = sqrt(p);
v1  = vr1 * r1   +   vt1 * cross(ih,r1);

vt2=vt1/r2mod;
vr2=-vr1+(vt1-vt2)/tan(theta/2.0);
v2=vr2*r2/r2mod+vt2*cross(ih,r2/r2mod);
if(e>1.0)
    a=-a;
end

%Subfunction that evaluates the time of flight as a function of x
function t=x2tof(x,s,c,lw,N)  

am=s/2.0;
a=am/(1-x*x);
if x<1.0 %ELLISSE
    beta=2.0*asin(sqrt((s-c)/2.0/a));
    if lw
        beta=-beta;
    end
    alfa=2.0*acos(x);
else   %IPERBOLE
    alfa=2.0*acosh(x);
    beta=2.0*asinh(sqrt((s-c)/(-2.0*a)));
    if lw
        beta=-beta;
    end
end
t=tofabn(a,alfa,beta,N);
%subfunction that evaluates the time of flight via Lagrange expression
function t=tofabn(sigma,alfa,beta,N)

if sigma>0.0
    t=sigma*sqrt(sigma)*((alfa-sin(alfa))-(beta-sin(beta))+N*2.0*pi);
else
    t=-sigma*sqrt(-sigma)*((sinh(alfa)-alfa)-(sinh(beta)-beta));
end
