%% Task 3.1b - Wanted output method

clear;
close all;
clc;

rng(100);

%% Generating samples

N = 500;

M1 = [5; 4];
M2 = [-5; 4];
M3 = [0; -5];

S1 = eye(2, 2);
S2 = eye(2, 2);
S3 = eye(2, 2);

X1 = mvnrnd(M1, S1, N)';
X2 = mvnrnd(M2, S2, N)';
X3 = mvnrnd(M3, S3, N)';

figure
hold all
scatter(X1(1, :), X1(2, :), 'r*')
scatter(X2(1, :), X2(2, :), 'g*')
scatter(X3(1, :), X3(2, :), 'b*')
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Location', 'southeast')
title('Odbirci klasa')
grid on
hold off

%% Wanted output method

gain1 = 1;
gain2 = 1;
gain3 = 1;

[V_12, v0_12] = wanted_output(X1, X2, gain1, gain2);
[V_23, v0_23] = wanted_output(X2, X3, gain2, gain3);
[V_31, v0_31] = wanted_output(X1, X3, gain1, gain3);

x = linspace(-8, 8, 1000);
y1 = -(V_12(1)*x + v0_12)/V_12(2);
y2 = -(V_23(1)*x + v0_23)/V_23(2);
y3 = -(V_31(1)*x + v0_31)/V_31(2);

figure
hold all
scatter(X1(1, :), X1(2, :), 'r*')
scatter(X2(1, :), X2(2, :), 'g*')
scatter(X3(1, :), X3(2, :), 'b*')
plot(x, y1, 'r--');
plot(x, y2, 'g--');
plot(x, y3, 'b--');
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasifikator 12', 'Klasifikator 23', 'Klasifikator 31', 'Location', 'southeast')
title('Odbirci klasa sa klasifikatorima dobijenim metodom zeljenog izlaza')
axis([-8 8 -10 8])
grid on
hold off

%% Confusion matrix

a1 = (y1(100) - y1(50))/(x(100) - x(50));
a2 = (y2(100) - y2(50))/(x(100) - x(50));
a3 = (y3(100) - y3(50))/(x(100) - x(50));
b1 = -a1*x(50) + y1(50);
b2 = -a2*x(50) + y2(50);
b3 = -a3*x(50) + y3(50);

N_errors_21 = errors_above(a1, b1, X2);
N_errors_12 = errors_underneath(a1, b1, X1);
N_errors_23 = errors_above(a2, b2, X3);
N_errors_32 = errors_underneath(a2, b2, X2);
N_errors_13 = errors_above(a3, b3, X3);
N_errors_31 = errors_underneath(a3, b3, X1);

% vrste su estimirane klase, kolone su tacne klase
conf_matrix = zeros(3, 3);
conf_matrix(1, 1) = length(X1) - N_errors_12 - N_errors_13;
conf_matrix(2, 2) = length(X2) - N_errors_21 - N_errors_23;
conf_matrix(3, 3) = length(X3) - N_errors_31 - N_errors_32;
conf_matrix(1, 2) = N_errors_12;
conf_matrix(1, 3) = N_errors_13;
conf_matrix(2, 1) = N_errors_21;
conf_matrix(2, 3) = N_errors_23;
conf_matrix(3, 1) = N_errors_31;
conf_matrix(3, 2) = N_errors_32;

disp('Konfuziona matrica: ');
disp(conf_matrix);
