%% Task 3.2 - Quadratic classsifier

clear;
close all;
clc;

rng(100);

%% Generating samples

N = 500;
 
R1 = rand(1,N);
F1 = rand(1,N)*2*pi;
X1 = [R1.*cos(F1); R1.*sin(F1)];

R2 = 4 + rand(1,N);
F2 = rand(1,N)*2*pi;
X2 = [R2.*cos(F2); R2.*sin(F2)];

figure
hold all
scatter(X1(1, :), X1(2, :), 'r*')
scatter(X2(1, :), X2(2, :), 'b*')
title('Odbirci klasa')
grid on
hold off

%% Wanted output method

gain1 = 1;
gain2 = 1;

[Q, V, v0] = wanted_output_sq(X1, X2, gain1, gain2);

x1 = -4.5:0.1:4.5;
x2 = -4.5:0.1:4.5;
h = zeros(length(x1),length(x2));

for i = 1:length(x1)
    for j = 1:length(x2)
        h(i,j) = v0+V(1)*x1(i)+V(2)*x2(j)+Q(1)*(x1(i))^2 + Q(2)*(x2(j))^2 + Q(3)*x1(i)*x2(j);
    end
end

figure
hold all
scatter(X1(1, :), X1(2, :), 'r*')
scatter(X2(1, :), X2(2, :), 'b*')
contour(x1, x2, h, [0 0], 'k', 'LineWidth',1);
legend('Klasa 1', 'Klasa 2', 'Klasifikator', 'Location', 'southeast')
title('Odbirci klasa sa kvadratnim klasifikatorom, dobijenim metodom zeljenog izlaza')
grid on
axis([-5 5 -5 5]);
hold off
