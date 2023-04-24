function [Y1_new, Y2_new, Y3_new, l] = ml_3(Y1, Y2, Y3)

T = 0.01;

% Inicijalna procena parametara

N1 = length(Y1(1,:));
N2 = length(Y2(1,:));
N3 = length(Y3(1,:));
N = N1 + N2 + N3;

P1 = N1/N;
P2 = N2/N;
P3 = N3/N;

M1 = mean(Y1,2);
M2 = mean(Y2,2);    
M3 = mean(Y3,2);

S1 = cov(Y1');
S2 = cov(Y2');
S3 = cov(Y3');


f1 = zeros(1, N);
f2 = zeros(1, N);
f3 = zeros(1, N);
f = zeros(1, N);

q1 = zeros(1, N);
q2 = zeros(1, N);
q3 = zeros(1, N);

q1_new = zeros(1, N);
q2_new = zeros(1, N);
q3_new = zeros(1, N);

Y = [Y1 Y2 Y3];

for i = 1:N
    f1(i) = 1/(2*pi*det(S1)^0.5)*exp(-0.5*(Y(:,i)-M1)'*inv(S1)*(Y(:,i)-M1));
    f2(i) = 1/(2*pi*det(S2)^0.5)*exp(-0.5*(Y(:,i)-M2)'*inv(S2)*(Y(:,i)-M2));
    f3(i) = 1/(2*pi*det(S3)^0.5)*exp(-0.5*(Y(:,i)-M3)'*inv(S3)*(Y(:,i)-M3));
    f(i) = P1*f1(i)+P2*f2(i)+P3*f3(i);
    
    q1(i) = P1*f1(i)/f(i);
    q2(i) = P2*f2(i)/f(i);
    q3(i) = P3*f3(i)/f(i);
end

end_alg = 0;
l=0;
while (end_alg == 0)
    
    l = l + 1;
    
    P1 = sum(q1)/length(q1);
    P2 = sum(q2)/length(q2);
    P3 = sum(q3)/length(q3);
    
    N1 = P1*N; 
    N2 = P2*N;
    N3 = P3*N;
    
    M1 = [0; 0]; M2 = [0; 0]; M3 = [0; 0];
    for i = 1:N
        M1 = M1 + q1(i)*Y(:,i);
        M2 = M2 + q2(i)*Y(:,i);
        M3 = M3 + q3(i)*Y(:,i);
    end
    M1 = M1/N1;
    M2 = M2/N2;
    M3 = M3/N3;
    
    S1 = zeros(2, 2);
    S2 = zeros(2, 2);
    S3 = zeros(2, 2);
    for i = 1:N
       S1 = S1 + q1(i)*(Y(:, i) - M1)*(Y(:, i) - M1)'; 
       S2 = S2 + q2(i)*(Y(:, i) - M2)*(Y(:, i) - M2)'; 
       S3 = S3 + q3(i)*(Y(:, i) - M3)*(Y(:, i) - M3)';
    end
    S1 = S1/N1;
    S2 = S2/N2;
    S3 = S3/N3;

    for i = 1:N
        f1(i) = 1/(2*pi*det(S1)^0.5)*exp(-0.5*(Y(:,i)-M1)'*inv(S1)*(Y(:,i)-M1));
        f2(i) = 1/(2*pi*det(S2)^0.5)*exp(-0.5*(Y(:,i)-M2)'*inv(S2)*(Y(:,i)-M2));
        f3(i) = 1/(2*pi*det(S3)^0.5)*exp(-0.5*(Y(:,i)-M3)'*inv(S3)*(Y(:,i)-M3));
        f(i) = P1*f1(i)+P2*f2(i)+P3*f3(i);

        q1_new(i) = P1*f1(i)/f(i);
        q2_new(i) = P2*f2(i)/f(i);
        q3_new(i) = P3*f3(i)/f(i);
    end

    dist1 = max(abs(q1-q1_new));
    dist2 = max(abs(q2-q2_new));
    dist3 = max(abs(q3-q3_new));
    
    d = max([dist1 dist2 dist3]);
    
    if d < T
        end_alg = 1;
    else
        clear q1 q2 q3
        q1 = q1_new;
        q2 = q2_new;
        q3 = q3_new;
    end
    
end

% Klasterizacija podataka

Y1_new = [];
Y2_new = [];
Y3_new = [];

for i = 1:N
    q = sort([q1_new(i), q2_new(i), q3_new(i)],'descend');
    if q(1) == q1_new(i)
        Y1_new = [Y1_new Y(:,i)];
    elseif q(1) == q2_new(i)
        Y2_new = [Y2_new Y(:,i)];
    else
        Y3_new = [Y3_new Y(:,i)];
    end
end

% Prikaz klastera

text = ['Potreban broj itercaija: ', num2str(l)];
disp(text);

figure
hold all
if ~isempty(Y1_new)
    plot(Y1_new(1, :), Y1_new(2, :), 'r*')
end
if ~isempty(Y2_new)
plot(Y2_new(1, :), Y2_new(2, :), 'b*')
end
if ~isempty(Y3_new)
plot(Y3_new(1, :), Y3_new(2, :), 'g*')
end
title('Odbirci klasa posle ml klasterizacije')
grid on
hold off

end

