%% Task 4.3 - Quadratic decomposition method
clear;
close all;
clc;

%% Generating samples

N = 500;
 
R1 = rand(1,N);
F1 = rand(1,N)*2*pi;
Y1 = [R1.*cos(F1); R1.*sin(F1)];

R2 = 3 + rand(1,N);
F2 = rand(1,N)*2*pi;
Y2 = [R2.*cos(F2); R2.*sin(F2)];

figure
hold all
scatter(Y1(1, :), Y1(2, :), 'r*')
scatter(Y2(1, :), Y2(2, :), 'b*')
title('Odbirci klasa')
grid on
hold off

%% Initial clustering

pom = rand(1, 2*N);
Y = [Y1 Y2];
Y1_start = [];
Y2_start = [];

for i=1:2*N
   if pom(i)<0.5
       Y1_start = [Y1_start Y(:, i)];
   else
       Y2_start = [Y2_start Y(:, i)];
   end
end

figure
hold all
scatter(Y1_start(1, :), Y1_start(2, :), 'r*')
scatter(Y2_start(1, :), Y2_start(2, :), 'b*')
legend('Klasa 1','Klasa 2');
title('Odbirci klasa posle pocetne klasterizacije')
grid on
hold off

%% Quadratic decomposition method

[Y1_decomp, Y2_decomp, l] = sq_decomp(Y1_start, Y2_start);

%% Dependence of the algorithm on initial clustering

pom = rand(1, 2*N);
Y = [Y1 Y2];
Y1_start1 = [];
Y2_start1 = [];
tr = linspace(0.2, 0.8, 10);
l = [];

for j=1:length(tr)
    for i=1:2*N
       if pom(i)<tr(j)
           Y1_start1 = [Y1_start1 Y(:, i)];
       else
           Y2_start1 = [Y2_start1 Y(:, i)];
       end
    end
    [Y1_decomp_i, Y2_decomp_i, l_i] = sq_decomp(Y1_start1, Y2_start1);
    l = [l l_i];
end

%% Mean number of iterations

l_mean = mean(l);
text = ['Srednji broj potrebnih iteracija: ', num2str(l_mean)];
disp(text);
