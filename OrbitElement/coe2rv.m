%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils                              
% Orbit elements ---> Cartesian position and velocity   
%                                                       
% Input: coe(6), mu, tol                                
%   coe(1): a, semi-major axis                          
%   coe(2): e, eccentricity                             
%   coe(3): i, orbit inclination                        
%   coe(4): Omega, right ascension of ascending node    
%   coe(5): omega, argument of pericentre               
%   coe(6): f, true anomaly                             
%   mu: gravitational coefficient                       
%   tol: tolerance, 1e-12 by default                    
%                                                       
% Output: CartesianR, CartesianV                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CartesianR, CartesianV] = coe2rv(coe, mu, tol)
% Number of inputs
if nargin < 3
    tol = 1e-12;
    if nargin < 2
        error("Not enough inputs.");
    end
end

% Orbit elements
% 6-dimensional
if length(coe) ~= 6
    error("Wrong length of orbit element vector. " + ...
        "The length should be 6.");
end

% Semi-major axis (a) & eccentricity (e)
% a>0, e>0
if coe(1) < 0 || coe(2) < 0
    error("Wrong semi-major axis or wrong eccentricity. " + ...
        "Both of them must be positive.");
end

% Orbit inclination
% i in [0,pi]
if coe(3) < 0 || coe(3) > pi
    error("Wrong orbit inclination. " + ...
        "Range: [0,pi]");
end

% To make sure r>0, we have "1+e*cos(f)>0"
% To avoid singular problem, add tolerance
if (coe(2) * cos(coe(6))) < -1.0 + tol
    error("Unsuitable input of eccentricity and true anomaly.");    
end

% Mu - must be positive
if mu < 0
    error("Wrong input mu. mu >= 0.");
end

% Orbit element list
a = coe(1);
e = coe(2);
i = coe(3);
Omega = coe(4);
omega = coe(5);
f = coe(6);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
% Semi-latus rectum
p = abs(a * (1 - e^2));

% Rotation matrix: 3-1-3 rotation
A3O = [cos(Omega), -sin(Omega), 0; 
       sin(Omega), cos(Omega),  0; 
       0,          0,           1];
A1i = [1, 0,      0;
       0, cos(i), -sin(i);
       0, sin(i), cos(i)];
A3o = [cos(omega), -sin(omega), 0;
       sin(omega), cos(omega),  0;
       0,          0,           1];
A = A3O * A1i * A3o;

% The direction of e-p coordinate
OX = [1; 0; 0];
ie = A * OX;
hd = [sin(i) * sin(Omega);
     -sin(i) * cos(Omega);
      cos(i)];
ip = cross(hd, ie);

% Output R and V (Cartesian coordinate)
r = p / (1 + e * cos(f));
CartesianR = r * cos(f) * ie + r * sin(f) * ip;
CartesianV = -sqrt(mu / p) * (sin(f) * ie - (e + cos(f)) * ip);
end