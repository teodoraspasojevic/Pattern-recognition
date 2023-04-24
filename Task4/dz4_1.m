%% Task 4.1 - C mean method

clear;
close all;
clc;

rng(100);

%% Generating samples

N = 500;
M1 = [3.5; 3.5];
M2 = [-3.5; 3.5];
M3 = [-3.5; -3.5];
M4 = [3.5; -3.5];
S1 = [1 0; 0 1];
S2 = [1 0; 0 1.5];
S3 = [1.2 0; 0 0.5];
S4 = [1 0; 0 1];

Y1 = mvnrnd(M1, S1, N)';
Y2 = mvnrnd(M2, S2, N)';
Y3 = mvnrnd(M3, S3, N)';
Y4 = mvnrnd(M4, S4, N)';

figure
hold all
scatter(Y1(1, :), Y1(2, :), 'r*')
scatter(Y2(1, :), Y2(2, :), 'b*')
scatter(Y3(1, :), Y3(2, :), 'g*')
scatter(Y4(1, :), Y4(2, :), 'y*')
title('Odbirci klasa')
grid on
hold off

%% Initial clustering

pom = rand(1, 4*N);
Y = [Y1 Y2 Y3 Y4];
Y1_start = [];
Y2_start = [];
Y3_start = [];
Y4_start = [];

for i=1:4*N
   if pom(i)<0.25
       Y1_start = [Y1_start Y(:, i)];
   elseif pom(i)<0.5
       Y2_start = [Y2_start Y(:, i)];
   elseif pom(i)<0.75
       Y3_start = [Y3_start Y(:, i)];
   else
       Y4_start = [Y4_start Y(:, i)];
   end
end

figure
hold all
scatter(Y1_start(1, :), Y1_start(2, :), 'r*')
scatter(Y2_start(1, :), Y2_start(2, :), 'b*')
scatter(Y3_start(1, :), Y3_start(2, :), 'g*')
scatter(Y4_start(1, :), Y4_start(2, :), 'y*')
legend('Klasa 1','Klasa 2','Klasa 3','Klasa 4');
title('Odbirci klasa posle pocetne klasterizacije')
grid on
hold off

%% C mean clustering

[Y1_cmean, Y2_cmean, Y3_cmean, Y4_cmean, l] = c_mean(Y1_start, Y2_start, Y3_start, Y4_start);

%% Dependence of the algorithm on initial clustering

pom = rand(1, 4*N);
Y = [Y1 Y2 Y3 Y4];
Y1_start = [];
Y2_start = [];
Y3_start = [];
Y4_start = [];
tr1 = linspace(0.2, 0.4, 10);
tr2 = linspace(0.4, 0.6, 10);
tr3 = linspace(0.6, 0.8, 10);
l = [];

for j=1:length(tr1)
    for i=1:4*N
       if pom(i)<tr1(j)
           Y1_start = [Y1_start Y(:, i)];
       elseif pom(i)<tr2(j)
           Y2_start = [Y2_start Y(:, i)];
       elseif pom(i)<tr3(j)
           Y3_start = [Y3_start Y(:, i)];
       else
           Y4_start = [Y4_start Y(:, i)];
       end
    end
    [Y1_cmean_i, Y2_cmean_i, Y3_cmean_i, Y4_cmean_i, l_i] = c_mean(Y1_start, Y2_start, Y3_start, Y4_start);
    l = [l l_i];
end

%% Mean number of iterations

l_mean = mean(l);
text = ['Srednji broj potrebnih iteracija: ', num2str(l_mean)];
disp(text);

%% Dependence of the algorithm on initial number of clusters

pom = rand(1, 4*N);
Y = [Y1 Y2 Y3 Y4];
Y1_start = [];
Y2_start = [];
Y3_start = [];

for i=1:4*N
   if pom(i)<0.33
       Y1_start = [Y1_start Y(:, i)];
   elseif pom(i)<0.66
       Y2_start = [Y2_start Y(:, i)];
   else
       Y3_start = [Y3_start Y(:, i)];
   end
end

figure
hold all
scatter(Y1_start(1, :), Y1_start(2, :), 'r*')
scatter(Y2_start(1, :), Y2_start(2, :), 'b*')
scatter(Y3_start(1, :), Y3_start(2, :), 'g*')
legend('Klasa 1','Klasa 2','Klasa 3');
title('Odbirci klasa posle pocetne klasterizacije')
grid on
hold off

[Y1_cmean3, Y2_cmean3, Y3_cmean3, l] = c_mean_3(Y1_start, Y2_start, Y3_start);

pom = rand(1, 4*N);
Y = [Y1 Y2 Y3 Y4];
Y1_start = [];
Y2_start = [];
Y3_start = [];
tr1 = linspace(0.2, 0.4, 10);
tr2 = linspace(0.4, 0.6, 10);
l = [];

for j=1:length(tr1)
    for i=1:4*N
       if pom(i)<tr1(j)
           Y1_start = [Y1_start Y(:, i)];
       elseif pom(i)<tr2(j)
           Y2_start = [Y2_start Y(:, i)];
       else
           Y3_start = [Y3_start Y(:, i)];
       end
    end
    [Y1_cmean_3_i, Y2_cmean_3_i, Y3_cmean_3_i, l_i] = c_mean_3(Y1_start, Y2_start, Y3_start);
    l = [l l_i];
end

l_mean = mean(l);
text = ['Srednji broj potrebnih iteracija: ', num2str(l_mean)];
disp(text);
