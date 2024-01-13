%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison between two transfers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
muSun = 1.327e11;
rEarth = 1.496e8;
rTerm = 2.279e8;

n = 100;
dvdList = zeros(n, 0);
dvhList = zeros(n, 0);
iList = linspace(0.2*pi, pi, n);
iEarth = 0;
for i=1:length(iList)
    iTerm = iList(i);

    [~, dvdList(i), ~, ~, ~, ~] = double_ellipse(rEarth, rTerm, ...
                                                iEarth, iTerm);
    [~, ~, dvhList(i), ~, ~, ~, ~] = hohmann(rEarth, rTerm, ...
                                                iEarth, iTerm);
end
f=figure;
plot(iList, dvhList, 'LineWidth', 1.5);hold on
plot(iList, dvdList, 'LineWidth', 1.5);
legend('Hohmann', 'Double-ellipse');
xlabel('i (rad)');
ylabel('dv (km/s)');
title('Comparison between Hohmann transfer and double-ellipse transfer');
saveas(f, 'comparison-i', 'fig');

%%
dvdList = zeros(n, 0);
dvhList = zeros(n, 0);
rList = linspace(1.5, 20, n);
iEarth = 0;
iTerm = 0;
for i=1:length(rList)
    rTerm = rList(i) * rEarth;

    [~, dvdList(i), ~, ~, ~, ~] = double_ellipse(rEarth, rTerm, ...
                                                iEarth, iTerm);
    [~, ~, dvhList(i), ~, ~, ~, ~] = hohmann(rEarth, rTerm, ...
                                                iEarth, iTerm);
end
f=figure;
plot(rList, dvhList, 'LineWidth', 1.5);hold on
plot(rList, dvdList, 'LineWidth', 1.5);
legend('Hohmann', 'Double-ellipse');
xlabel('r (magnification)');
ylabel('dv (km/s)');
title('Comparison between Hohmann transfer and double-ellipse transfer');
saveas(f, 'comparison-r', 'fig');



