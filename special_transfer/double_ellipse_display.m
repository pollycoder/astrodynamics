%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Result display for Double-ellipse transfer
% Input:
%   The output of double-ellipse solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function double_ellipse_display(output, title)
fprintf(title);
fprintf("\n");
fprintf("dv1 = %f\n", output(1));
fprintf("dv2 = %f\n", output(2));
fprintf("dv3 = %f\n", output(3));
fprintf("dv = %f\n", output(4));
fprintf("dt = %f\n", output(5));
fprintf("a1 = %f\n", output(6));
fprintf("e1 = %f\n", output(7));
fprintf("i1 = %f\n", output(8));
fprintf("a2 = %f\n", output(9));
fprintf("e2 = %f\n", output(10));
fprintf("i2 = %f\n", output(11));
end