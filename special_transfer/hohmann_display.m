%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Result display for Hohmann transfer
% Input:
%   The output of Hohmann solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hohmann_display(output, title)
fprintf(title);
fprintf("\n");
fprintf("dv1 = %f\n", output(1));
fprintf("dv2 = %f\n", output(2));
fprintf("dv = %f\n", output(3));
fprintf("dt = %f\n", output(4));
fprintf("a = %f\n", output(5));
fprintf("e = %f\n", output(6));
fprintf("i = %f\n", output(7));
end