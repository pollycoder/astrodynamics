function [e1,e2,e3]=transformation(x,y,z)
%将地固直角坐标系中的直角坐标为(x, y, z)转换到球坐标系中
e1=[x,y,z]'/sqrt(x^2+y^2+z^2);
e2=[-y,x,0]'/sqrt(x^2+y^2);
e3=cross(e1,e2);
end