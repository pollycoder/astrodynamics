%%%%%%%%%%%%%%%%%%%%%%%%
% Object function
%%%%%%%%%%%%%%%%%%%%%%%%

% X(0):start time
% X(1):end time
function result=impulse_obj(X)
day = 86400;
muSun = 1.327e11;
aEarth = 1.496e8;
aMars = 2.279e8;
iEarth = 0;
iMars = deg2rad(1.85);
eEarth = 0.0167;
eMars = 0.0549;

% Set the time
delta_X=X(2)-X(1);                      % Normalized transfer_time
if delta_X<=0 || delta_X>5              % Penalty function
    result=1E20;
    return
end
wait_time=X(1)*day*365;                 % Waiting time (s)
transfer_time=delta_X*day*365;          % Transfer time (s)
total_time=X(2)*day*365;                % Total time (s)

f0 = 0;
ft = 2/ 3 * pi;
f0t = f0dt2ft(f0, wait_time, aEarth, 0, muSun);
coeEarth = [aEarth, eEarth, iEarth, 0, 0, f0t];
ftt = f0dt2ft(ft, total_time, aMars, 0, muSun);
coeMars = [aMars, eMars, iMars,0, 0, ftt];

% Initial position
[r0,v0]=coe2rv(coeEarth, muSun);

% Final position
[rf,vf]=coe2rv(coeMars, muSun);

% Lambert
[v1,v2,~,~,~,~]=LambSol(r0,rf,transfer_time,muSun);

dv1=delta_v(v1,v0);
dv2=delta_v(v2,vf);
result=dv1+dv2;
end

function dv=delta_v(v1,v2)
v=v1-v2;
dv=norm(v,2);
end