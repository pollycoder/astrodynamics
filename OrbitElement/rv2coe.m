%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbit Calculation Utils                              
% Cartesian position and velocity ---> Orbit elements  
%                                                       
% Input: CartesianR, CartesianV, mu, tol                                
%   CartesianR: position in Cartesian coordinate system
%   CartesianV: velocity in Cartesian coordiante system
%   mu: gravitational coefficient                       
%   tol: tolerance, 1e-12 by default                    
%                                                       
% Output: coe(6)  
%   coe(1): a, semi-major axis                          
%   coe(2): e, eccentricity                             
%   coe(3): i, orbit inclination                        
%   coe(4): Omega, right ascension of ascending node    
%   coe(5): omega, argument of pericentre               
%   coe(6): f, true anomaly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coe = rv2coe(CartesianR, CartesianV, mu)
if nargin == 2
    mu = 398600;
end

r = norm(CartesianR);
v = norm(CartesianV);

OX = [1; 0; 0];
OY = [0; 1; 0];
OZ = [0; 0; 1];

% Semi-major axis (a)
Em = 0.5 * v^2 - mu / r;
a = abs(-mu / (2 * Em));
coe(1) = a;

% Eccentric vector
h = cross(CartesianR, CartesianV);
hNorm = norm(h);
e = 1 / mu * cross(CartesianV, h) - CartesianR / r;
eNorm = norm(e);
coe(2) = eNorm;

% Orbit inclination
cosi = sum(h .* OZ) / hNorm;
i = acos(cosi);
coe(3) = i;

% Right ascension of ascending node
hVector = h / hNorm;
ON = cross(OZ, hVector);                        % Ascending pitch line
if norm(ON) == 0
    ON = OX;
end
ON = ON / norm(ON);


Omega0 = acos(ON' * OX);
temp = ON' * OY;                                % Judge whether Omega>pi                               
if temp > 0                                     % ON*OY>0, Omega<pi       
    Omega = Omega0;
else                                            % ON*OY<0, Omega>pi
    Omega = 2 * pi - Omega0;
end
coe(4) = Omega;

% Argument of pericentre
eVector = e / eNorm;
temp = h' * cross(ON, eVector);                 % Judge whether omega>pi
omega0 = acos(ON' * eVector);
if temp > 0                                     % omega<pi
    omega = omega0;
else                                            % omega>pi
    omega = 2 * pi - omega0;
end
coe(5) = omega;

% True anomaly
temp = CartesianR' * CartesianV;                % Judge whether f>pi
f0 = acos(eVector' * CartesianR / r);
if temp > 0                                     % f<pi
    f = f0;
else
    if eNorm >= 1                               % Hyperbola and parabola
        f = -f0;
    else                                        % Ellipse and circle
        f = 2 * pi - f0;  
    end
end
coe(6) = f;

end