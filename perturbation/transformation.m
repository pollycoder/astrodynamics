function [e1,e2,e3]=transformation(x,y,z)
%���ع�ֱ������ϵ�е�ֱ������Ϊ(x, y, z)ת����������ϵ��
e1=[x,y,z]'/sqrt(x^2+y^2+z^2);
e2=[-y,x,0]'/sqrt(x^2+y^2);
e3=cross(e1,e2);
end