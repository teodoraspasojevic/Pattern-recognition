%% Task 4.2 - Maximum likelihood algorithm

clear;
close all;
clc;

rng(100);

%% Generating samples

N = 500;
M1 = [1 1]';
M2 = [1 5]';
M3 = [5 5]';
M4 = [5 1]';
S1 = eye(2, 2)*0.5;
S2 = eye(2, 2)*0.5;
S3 = eye(2, 2)*0.5;
S4 = eye(2, 2)*0.5;

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
title('Odbirci klasa posle pocetne klasterizacije')
grid on
hold off

%% ML clustering

[Y1_new, Y2_new, Y3_new, Y4_new] = ml(Y1_start, Y2_start, Y3_start, Y4_start);

%% Adding prior knowledge to inital clustering

Y1_start1 = Y1(:, 1:N/4);
Y2_start1 = Y2(:, 1:N/4);
Y3_start1 = Y3(:, 1:N/4);
Y4_start1 = Y4(:, 1:N/4);

for i=N/4+1:N
   if pom(i)<0.25
       Y1_start1 = [Y1_start1 Y1(:, i)];
   elseif pom(i)<0.5
       Y2_start1 = [Y2_start1 Y1(:, i)];
   elseif pom(i)<0.75
       Y3_start1 = [Y3_start1 Y1(:, i)];
   else
       Y4_start1 = [Y4_start1 Y1(:, i)];
   end
end
for i=N/4+1:N
   if pom(i + N)<0.25
       Y1_start1 = [Y1_start1 Y2(:, i)];
   elseif pom(i)<0.5
       Y2_start1 = [Y2_start1 Y2(:, i)];
   elseif pom(i)<0.75
       Y3_start1 = [Y3_start1 Y2(:, i)];
   else
       Y4_start1 = [Y4_start1 Y2(:, i)];
   end
end
for i=N/4+1:N
   if pom(i + 2*N)<0.25
       Y1_start1 = [Y1_start1 Y3(:, i)];
   elseif pom(i)<0.5
       Y2_start1 = [Y2_start1 Y3(:, i)];
   elseif pom(i)<0.75
       Y3_start1 = [Y3_start1 Y3(:, i)];
   else
       Y4_start1 = [Y4_start1 Y3(:, i)];
   end
end
for i=N/4+1:N
   if pom(i + 3*N)<0.25
       Y1_start1 = [Y1_start1 Y4(:, i)];
   elseif pom(i)<0.5
       Y2_start1 = [Y2_start1 Y4(:, i)];
   elseif pom(i)<0.75
       Y3_start1 = [Y3_start1 Y4(:, i)];
   else
       Y4_start1 = [Y4_start1 Y4(:, i)];
   end
end

[Y1_new1, Y2_new1, Y3_new1, Y4_new1] = ml(Y1_start1, Y2_start1, Y3_start1, Y4_start1);

%% Dependence of the algorithm on initial clustering

pom = rand(1, 4*N);
Y = [Y1 Y2 Y3 Y4];
Y1_start2 = [];
Y2_start2 = [];
Y3_start2 = [];
Y4_start2 = [];
tr1 = linspace(0.2, 0.4, 10);
tr2 = linspace(0.4, 0.6, 10);
tr3 = linspace(0.6, 0.8, 10);
l = [];

for j=1:length(tr1)
    
    Y1_start2 = Y1(:, 1:N/4);
    Y2_start2 = Y2(:, 1:N/4);
    Y3_start2 = Y3(:, 1:N/4);
    Y4_start2 = Y4(:, 1:N/4);

    for i=N/4+1:N
        if pom(i)<tr1(j)
            Y1_start2 = [Y1_start2 Y1(:, i)];
        elseif pom(i)<tr2(j)
            Y2_start2 = [Y2_start2 Y1(:, i)];
        elseif pom(i)<tr3(j)
            Y3_start2 = [Y3_start2 Y1(:, i)];
        else
            Y4_start2 = [Y4_start2 Y1(:, i)];
        end
    end
    for i=N/4+1:N
        if pom(i + N)<tr1(j)
            Y1_start2 = [Y1_start2 Y2(:, i)];
        elseif pom(i)<tr2(j)
            Y2_start2 = [Y2_start2 Y2(:, i)];
        elseif pom(i)<tr3(j)
            Y3_start2 = [Y3_start2 Y2(:, i)];
        else
            Y4_start2 = [Y4_start2 Y2(:, i)];
        end
    end
    for i=N/4+1:N
        if pom(i + 2*N)<tr1(j)
            Y1_start2 = [Y1_start2 Y3(:, i)];
        elseif pom(i)<tr2(j)
            Y2_start2 = [Y2_start2 Y3(:, i)];
        elseif pom(i)<tr3(j)
            Y3_start2 = [Y3_start2 Y3(:, i)];
        else
            Y4_start2 = [Y4_start2 Y3(:, i)];
        end
    end
    for i=N/4+1:N
        if pom(i + 3*N)<tr1(j)
            Y1_start2 = [Y1_start2 Y4(:, i)];
        elseif pom(i)<tr2(j)
            Y2_start2 = [Y2_start2 Y4(:, i)];
        elseif pom(i)<tr3(j)
            Y3_start2 = [Y3_start2 Y4(:, i)];
        else
            Y4_start2 = [Y4_start2 Y4(:, i)];
        end
    end
    [Y1_new_i, Y2_new_i, Y3_new_i, Y4_new_i, l_i] = ml(Y1_start2, Y2_start2, Y3_start2, Y4_start2);
    l = [l l_i];
end

%% Mean number of iterations

l_mean = mean(l);
text = ['Srednji broj potrebnih iteracija: ', num2str(l_mean)];
disp(text);

%% Dependence of the algorithm on initial number of clusters

pom = rand(1, 4*N);
Y = [Y1 Y2 Y3 Y4];
Y1_start3 = [];
Y2_start3 = [];
Y3_start3 = [];

for i=1:4*N
   if pom(i)<0.33
       Y1_start3 = [Y1_start3 Y(:, i)];
   elseif pom(i)<0.66
       Y2_start3 = [Y2_start3 Y(:, i)];
   else
       Y3_start3 = [Y3_start3 Y(:, i)];
   end
end

figure
hold all
scatter(Y1_start3(1, :), Y1_start3(2, :), 'r*')
scatter(Y2_start3(1, :), Y2_start3(2, :), 'b*')
scatter(Y3_start3(1, :), Y3_start3(2, :), 'g*')
legend('Klasa 1','Klasa 2','Klasa 3');
title('Odbirci klasa posle pocetne klasterizacije')
grid on
hold off

[Y1_new3, Y2_new3, Y3_new3, l] = ml_3(Y1_start3, Y2_start3, Y3_start3);
 
pom = rand(1, 4*N);
Y = [Y1 Y2 Y3 Y4];
Y1_start4 = [];
Y2_start4 = [];
Y3_start4 = [];
tr1 = linspace(0.2, 0.4, 10);
tr2 = linspace(0.4, 0.6, 10);
l = [];

for j=1:length(tr1)
    
    Y1_start4 = Y1(:, 1:N/4);
    Y2_start4 = Y2(:, 1:N/4);
    Y3_start4 = Y3(:, 1:N/4);

    for i=N/4+1:N
        if pom(i)<tr1(j)
            Y1_start4 = [Y1_start4 Y1(:, i)];
        elseif pom(i)<tr2(j)
            Y2_start4 = [Y2_start4 Y1(:, i)];
        else
            Y3_start4 = [Y3_start4 Y1(:, i)];
        end
    end
    for i=N/4+1:N
        if pom(i + N)<tr1(j)
            Y1_start4 = [Y1_start4 Y2(:, i)];
        elseif pom(i)<tr2(j)
            Y2_start4 = [Y2_start4 Y2(:, i)];
        else
            Y3_start4 = [Y3_start4 Y2(:, i)];
        end
    end
    for i=N/4+1:N
        if pom(i + 2*N)<tr1(j)
            Y1_start4 = [Y1_start4 Y3(:, i)];
        elseif pom(i)<tr2(j)
            Y2_start4 = [Y2_start4 Y3(:, i)];
        else
            Y3_start4 = [Y3_start4 Y3(:, i)];
        end
    end
    for i=N/4+1:N
        if pom(i + 3*N)<tr1(j)
            Y1_start4 = [Y1_start4 Y4(:, i)];
        elseif pom(i)<tr2(j)
            Y2_start4 = [Y2_start4 Y4(:, i)];
        else
            Y3_start4 = [Y3_start4 Y4(:, i)];
        end
    end
    [Y1_new_i, Y2_new_i, Y3_new_i, l_i] = ml_3(Y1_start2, Y2_start2, Y3_start2);
    l = [l l_i];
end

l_mean = mean(l);
text = ['Srednji broj potrebnih iteracija: ', num2str(l_mean)];
disp(text);
