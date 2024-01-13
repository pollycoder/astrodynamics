%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remember to run this section first before 
% you run any other sections !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
au = 1.496e8;
muSun = 1.327e11;
rEarth = 1 * au;
rTerm = 20 * rEarth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Double-ellipse Transfer (AU, s)
% Earth and Mars
% Coplanar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iEarth = 0;
iTerm = 0;

% Double-ellipse transfer
[DVd, dvd, dtd, atd, etd, itd] = double_ellipse(rEarth, rTerm, ...
                                                iEarth, iTerm, muSun);
output = [DVd(1), DVd(2), DVd(3), dvd, dtd, ...
          atd(1), atd(2), etd(1), etd(2), itd(1), itd(2)];
title = 'DE - coplanar';
double_ellipse_display(output, title);

% Hohmann transfer
[dv1h, dv2h, dvh, dth, ah, eh, ih] = hohmann(rEarth, rTerm, ...
                                             iEarth, iTerm, muSun);
output = [dv1h, dv2h, dvh, dth, ah, eh, ih];
title = 'Hohmann - coplanar';
hohmann_display(output, title);

% Comparison
if dvd < dvh
    fprintf("Double-ellipse transfer is better than Hohmann transfer.\n");
else
    fprintf("Hohmann transfer is good enough.\n");
end

figure
% Plot
% Hohmann
f1=figure(1);
coeD1 = [rEarth, iEarth];
coeA1 = [rTerm, iTerm];
coeH1 = [ah, eh, ih];
hohmann_plot(coeD1, coeA1, coeH1, muSun, 'coplanar');

% Double-ellipse
f2=figure(2);
coet1 = [atd(1), etd(1), itd(1)];
coet2 = [atd(2), etd(2), itd(2)];
double_ellipse_plot(coeD1, coeA1, coet1, coet2, muSun, 'coplanar');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Double-ellipse Transfer (AU, s)
% Earth and Mars
% Non-coplanar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rTerm = 2 * rEarth;
iEarth = 0;
iTerm = 0.25 * pi;

% Double-ellipse transfer
[DVd, dvd, dtd, atd, etd, itd] = double_ellipse(rEarth, rTerm, ...
                                                iEarth, iTerm, muSun);
output = [DVd(1), DVd(2), DVd(3), dvd, dtd, ...
          atd(1), atd(2), etd(1), etd(2), itd(1), itd(2)];
title = 'DE - noncoplanar';
double_ellipse_display(output, title);

% Hohmann transfer
[dv1h, dv2h, dvh, dth, ah, eh, ih] = hohmann(rEarth, rTerm, ...
                                             iEarth, iTerm, muSun);
output = [dv1h, dv2h, dvh, dth, ah, eh, ih];
title = 'Hohmann - noncoplanar';
hohmann_display(output, title);

% Comparison
if dvd < dvh
    fprintf("Double-ellipse transfer is better than Hohmann transfer.\n");
else
    fprintf("Hohmann transfer is good enough.\n");
end


% Plot
% Hohmann
f3=figure(3);
coeD1 = [rEarth, iEarth];
coeA1 = [rTerm, iTerm];
coeH1 = [ah, eh, ih];
hohmann_plot(coeD1, coeA1, coeH1, muSun, 'non-coplanar');

% Double-ellipse
f4=figure(4);
coet1 = [atd(1), etd(1), itd(1)];
coet2 = [atd(2), etd(2), itd(2)];
double_ellipse_plot(coeD1, coeA1, coet1, coet2, muSun, 'non-coplanar');

% Save file
saveas(f1, 'hohmann-coplanar', 'fig');
saveas(f2, 'double-ellipse-coplanar', 'fig');
saveas(f3, 'hohmann-noncoplanar', 'fig');
saveas(f4, 'double-ellipse-noncoplanar', 'fig');