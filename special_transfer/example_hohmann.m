%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remember to run this section first before 
% you run any other sections !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
muSun = 1.327e11;
rEarth = 1.496e8;
rMars = 4.37e8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hohmann trasfer (km, s)
% Earth and Mars
% Coplanar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iEarth = 0;
iMars = 0;
[dv1, dv2, dv, dt, ah, eh, ih] = hohmann(rEarth, rMars, iEarth, iMars);

figure
coeD1 = [rEarth, iEarth];
coeA1 = [rMars, iMars];
coeH1 = [ah, eh, ih];
output1 = [dv1, dv2, dv, dt, ah, eh, ih];
title = "coplanar";
hohmann_display(output1, title);
hohmann_plot(coeD1, coeA1, coeH1, muSun, title);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hohmann trasfer (km, s)
% Earth and Mars
% Non-coplanar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iMars = deg2rad(1.850);
[dv1, dv2, dv, dt, ah, eh, ih] = hohmann(rEarth, rMars, iEarth, iMars);

figure
coeD2 = [rEarth, iEarth];
coeA2 = [rMars, iMars];
coeH2 = [ah, eh, ih];
output2 = [dv1, dv2, dv, dt, ah, eh, ih];
title = "noncoplanar";
hohmann_display(output2, title);
hohmann_plot(coeD2, coeA2, coeH2, muSun, title);