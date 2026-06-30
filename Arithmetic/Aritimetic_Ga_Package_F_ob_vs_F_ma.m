%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements the technique proposed in the article titled:
% "A Novel Transformation Technique for Solving Highly Linear Systems of Equations via Evolutionary Algorithms".
% This program uses the GA-octave package to run genetic algorithms.
% To run this program correctly, ensure the GA-octave package is installed.
% and this file is saved in the correct folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all  % Clear all variables from the workspace
pkg load statistics

num_runs = 20;  % Number of times each sub-block of the algorithm will run

% Define Bounds for the First Block (10 variables)
lb1 = -ones(1, 10);  % Generates a row of ten -1s
ub1 = ones(1, 10);   % Generates a row of ten 1s

% Define Bounds for the Second Block (5 variables)
lb2 = -ones(1, 5);   % Generates a row of five -1s
ub2 = ones(1, 5);    % Generates a row of five 1s


% =================================================================
% DEFINITION OF OBJECTIVE FUNCTIONS
% =================================================================

function z = funt(x)
    y = zeros(10, 1);  % Initialize y with zeros
    y(1) = x(6) - 0.16275449 * (x(1) * x(5) * x(9)) - 0.37842197;
    y(2) = x(7) - 0.15585316 * (x(3) * x(1) * x(9)) - 0.19807914;
    y(3) = x(8) - 0.19950920 * (x(3) * x(9) * x(2)) - 0.44166728;
    y(4) = x(9) - 0.18922793 * (x(4) * x(8) * x(5)) - 0.14654113;
    y(5) = x(5) - 0.19612740 * (x(5) * x(9) * x(4)) - 0.34504906;
    y(6) = x(1) - 0.18324757 * (x(7) * x(2) * x(10)) - 0.25428722;
    y(7) = x(2) - 0.16955071 * (x(1) * x(6) * x(5)) - 0.27162577;
    y(8) = x(3) - 0.21180486 * (x(6) * x(8) * x(4)) - 0.42937161;
    y(9) = x(4) - 0.17081208 * (x(1) * x(3) * x(9)) - 0.07056438;
    y(10) = x(10) - 0.21466544 * (x(7) * x(4) * x(1)) - 0.42651102;
    z = max(abs(y));  % Objective is the maximum absolute value of y
endfunction

function z = fun(x)
    S = [0.37842197 0.19807914 0.44166728 0.14654113 0.34504906];
    SE = [0.25428722, 0.27162577, 0.42937161, 0.07056438, 0.42651102];

    zh = zeros(5, 5);  % Initialize zh matrix with zeros
    zh(1, 1) = 1; zh(1, 4) = -0.16275449 * (x(1) * x(5));
    zh(2, 2) = 1; zh(2, 4) = -0.15585316 * (x(3) * x(1));
    zh(3, 3) = 1; zh(3, 4) = -0.19950920 * (x(3) * x(2));
    zh(4, 3) = -0.18922793 * (x(4) * x(5)); zh(4, 4) = 1;
    zh(5, 4) = -0.19612740 * (x(4) * x(5)); zh(5, 5) = 1;

    w = zeros(1, 10);  % Initialize w with zeros
    w(6:10) = (zh \ S')';  % Solve for w

    y = zeros(5, 1);  % Initialize y with zeros
    y(1) = x(1) - 0.18324757 * (w(7) * x(2) * w(10));
    y(2) = x(2) - 0.16955071 * (x(1) * w(6) * x(5));
    y(3) = x(3) - 0.21180486 * (w(6) * w(8) * x(4));
    y(4) = x(4) - 0.17081208 * (x(1) * x(3) * w(9));
    y(5) = x(5) - 0.21466544 * (w(7) * x(4) * x(1));

    z = sum(abs(y' - SE));  % Objective is the sum of absolute differences
endfunction


% =================================================================
% ALGORITHM EXECUTION
% =================================================================

% --- First Run Block (F_ma) ---
tic  % Start timer
for i = 1: num_runs
    [root, fval] = ga(@funt, 10, [], [], [], [], lb1, ub1);
    ST(i, :) = [fval, root];  % Store results
endfor
time1 = toc;  % Stop timer
TMt1 = time1 / num_runs;  % Calculate average time per generation

% --- Second Run Block (F_ob) ---
tic  % Start timer
for i = 1:num_runs
    [root, fval] = ga(@fun, 5, [], [], [], [], lb2, ub2);
    SF(i, :) = [fval, root];  % Store results
endfor
time_2 = toc;  % Stop timer
TMt2 = time_2 / num_runs;  % Calculate average time per generation


% =================================================================
% RESULTS AND PLOTTING
% =================================================================

xp = 1 : num_runs;  % X-axis values
ypt = ST(:, 1);  % Y-axis values for first objective function (F_ma)
yp = SF(:, 1);   % Y-axis values for second objective function (F_ob)

plot(xp, yp, ':bs', 'LineWidth', 3, xp, ypt, ':ko', 'LineWidth', 2);
legend('Objective function F_{ob}' , 'Objective function F_{ma}');  % Add legend
xlabel('Solutions');  % X-axis label
ylabel('Value objective function');  % Y-axis label

% --- Statistical Non-Parametric Test ---
disp('--------------------------------------------------')
disp('Calculating the Wilcoxon Rank-Sum Test (Mann-Whitney)...')
p_value = ranksum(yp, ypt)
disp('--------------------------------------------------')
