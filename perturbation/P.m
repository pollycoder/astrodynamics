function Pnm=P(n,m,u)
%归一化的缔合勒让德多项式
if(n==0&&m==0)
    Pnm=1;
elseif(n==1&&m==0)
    Pnm=sqrt(3)*u;
elseif(n==1&&m==1)
    Pnm=sqrt(3)*sqrt(1-u^2);
elseif(n<m)
    Pnm=0;
elseif(n==m)
    Pnm=sqrt(3)*sqrt(1-u^2);
    for i=2:n
        Pnm=Pnm*sqrt((2*i+1)/(2*i))*sqrt(1-u^2);
    end
elseif(n>m)
    Pnm0=P(m,m,u);
    Pnm10=P(1,0,u);
    Pnmup = 0; Pnm =0;
    for i=m:n
        if i-2<m
            Pnm2=0;
        elseif i-2==m
            Pnm2=Pnm0*sqrt(((2*i+1)*(i-1+m)*(i-1-m))/((2*i-3)*(i+m)*(i-m)));
        elseif i-2==1&&m==0
            Pnm2=Pnm10*sqrt(((2*i+1)*(i-1+m)*(i-1-m))/((2*i-3)*(i+m)*(i-m)));
        else
            Pnm2=Pnmup*sqrt(((2*i+1)*(i-1+m)*(i-1-m))/((2*i-3)*(i+m)*(i-m)));
        end
        if i-1<m
            Pnm1=0;
        elseif i-1==m
            Pnm1=Pnm0*sqrt((2*i+1)*(2*i-1)/((i+m)*(i-m)))*u;
        elseif i-1==1&&m==0
            Pnm1=Pnm10*sqrt((2*i+1)*(2*i-1)/((i+m)*(i-m)))*u;
        else
            Pnm1=Pnm*sqrt((2*i+1)*(2*i-1)/((i+m)*(i-m)))*u;
        end
        Pnmup = Pnm;
        Pnm=Pnm1-Pnm2;
    end
end
end