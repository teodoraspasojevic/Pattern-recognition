%% Task 2 - Classification based on hypothesis testing

clear;
close all;
clc;

rng(100);

%% Generating samples

N = 500;

M11 = [4.5; 3]; M12 = [4.5; 5];
M21 = [1; 3]; M22 = [1; 5];

S11 = eye(2, 2)*0.5; S12 = eye(2, 2)*0.4;
S21 = eye(2, 2)*0.6; S22 = eye(2, 2)*0.3;

X11 = mvnrnd(M11, S11, N)';
X12 = mvnrnd(M12, S12, N)';
X21 = mvnrnd(M21, S21, N)';
X22 = mvnrnd(M22, S22, N)';

P11 = 0.7; P12 = 0.3;
P21 = 0.5; P22 = 0.5;

pom1 = rand(1, N);
pom2 = rand(1, N);
X1 = (pom1 < P11).*X11 + (pom1 >= P11).*X12;
X2 = (pom2 < P21).*X21 + (pom2 >= P21).*X22;

figure
hold all
scatter(X1(1, :), X1(2, :), 'r*')
scatter(X2(1, :), X2(2, :), 'b*')
legend('Klasa 1', 'Klasa 2', 'Location', 'southeast')
title('Odbirci klasa')
grid on
hold off

%% Probability density function and histogram of classes

x = linspace(-1, 7, 100);
y = linspace(-1, 7, 100);
f1 = zeros(length(x),length(y));
f2 = zeros(length(x),length(y));

const1 = 1/(2*pi*det(S11)^(0.5));
const2 = 1/(2*pi*det(S12)^(0.5));
const3 = 1/(2*pi*det(S21)^(0.5));
const4 = 1/(2*pi*det(S22)^(0.5));

for i = 1:length(x)
    for j = 1:length(y)
        X = [x(i) y(j)]';
        f11 = const1*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
        f12 = const2*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
        f21 = const3*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
        f22 = const4*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
        f1(i,j) = P11*f11 + P12*f12;
        f2(i,j) = P21*f21 + P22*f22;
   end
end

figure
hold all
mesh(x, y, f1');
mesh(x, y, f2');
legend('Klasa 1', 'Klasa 2')
title('Funkcije gustine verovatnoce klasa')
hold off

figure
histogram2(X1, X2, 40, 'Normalization', 'pdf')
title('Histogram klasa');

%% Bayesian classifier

P1 = 0.5; P2 = 0.5;

T = log(P1/P2);
h = -log(f1./f2);

figure
hold all
scatter(X1(1, :), X1(2, :), 'r*')
scatter(X2(1, :), X2(2, :), 'b*')
contour(x, y, h', [T, T], 'k', 'LineWidth', 1)
legend('Klasa 1', 'Klasa 2', 'Bajesov klasifikator', 'Location', 'southeast')
title('Bajesov klasifikator minimalne verovatnoce greske')

%% Computation of the error of classifier

% confusion matrix
conf_matr = zeros(2,2);

const1 = 1/(2*pi*det(S11)^(0.5));
const2 = 1/(2*pi*det(S12)^(0.5));
const3 = 1/(2*pi*det(S21)^(0.5));
const4 = 1/(2*pi*det(S22)^(0.5));

for i=1:N
    X = X1(:, i);
    f11 = const1*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
    f12 = const2*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
    f21 = const3*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
    f22 = const4*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
    f1_i = P11*f11 + P12*f12;
    f2_i = P21*f21 + P22*f22;
    h = -log(f1_i/f2_i);
    if h<0
        conf_matr(1,1) = conf_matr(1,1) + 1;
    else 
        conf_matr(2,1) = conf_matr(2,1) + 1;
    end
end

for i=1:N
    X = X2(:, i);
    f11 = const1*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
    f12 = const2*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
    f21 = const3*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
    f22 = const4*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
    f1_i = P11*f11 + P12*f12;
    f2_i = P21*f21 + P22*f22;
    h = -log(f1_i/f2_i);
    if h<0
        conf_matr(1,2) = conf_matr(1,2) + 1; % greska
    else 
        conf_matr(2,2) = conf_matr(2,2) + 1; % pogodak
    end
end

disp('Konfuziona matrica:');
disp(conf_matr);

error = (conf_matr(1,2)+conf_matr(2,1))/sum(sum(conf_matr));
error1 = conf_matr(1,2)/N; % verovatnoca greske drugog tipa
error2 = conf_matr(2,1)/N; % verovatnoca greske prvog tipa

text1 = ['Verovatnoca greske prvog tipa - iz konfuzione matrice: ', num2str(error1)];
text2 = ['Verovatnoca greske drugog tipa - iz konfuzione matrice: ', num2str(error2)];

disp(text1);
disp(text2);

% theoretical approach
x1 = -1:0.1:7;
y1 = -1:0.1:7;

const1 = 1/(2*pi*det(S11)^(0.5));
const2 = 1/(2*pi*det(S12)^(0.5));
const3 = 1/(2*pi*det(S21)^(0.5));
const4 = 1/(2*pi*det(S22)^(0.5));

err1 = 0;
err2 = 0;

for i = 1:length(x1)
    for j = 1:length(y1)
        X = [x1(i) y1(j)]';
        f11 = const1*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
        f12 = const2*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
        f21 = const3*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
        f22 = const4*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
        f1_i = P11*f11 + P12*f12;
        f2_i = P21*f21 + P22*f22;
        h = -log(f1_i/f2_i);
        if h < 0
            err2 = err2 + f2_i*0.1*0.1; % greska - integralimo f2 po oblasti L1
        else
            err1 = err1 + f1_i*0.1*0.1; % greska - integralimo f1 po oblasti L2
        end
   end
end

text1 = ['Verovatnoca greske prvog tipa - teorijski pristup: ', num2str(err1)];
text2 = ['Verovatnoca greske drugog tipa - teorijski pristup: ', num2str(err2)];

disp(text1);
disp(text2);

%% Minimum price classifier

c11 = 0; c22 = 0;
c12 = 1; c21 = 10;

T = log((c21 - c11)*P1/(c12 - c22)/P2);
h = -log(f1./f2);

figure
hold all
scatter(X1(1, :), X1(2, :), 'r*')
scatter(X2(1, :), X2(2, :), 'b*')
contour(x, y, h', [T, T], 'k', 'LineWidth', 1)
legend('Klasa 1', 'Klasa 2', 'Klasifikator min cene', 'Location', 'southeast')
title('Klasifikator minimalne cene')

%% Neyman-Pearson classifier

h = -log(f1./f2);

mi= 0.1:0.1:10;
error2_v = zeros(1, length(mi));
for i=1:length(mi)
    for j = 1:length(x)-1
        for k = 1:length(y)-1
            if (h(j, k) < -log(mi(i)))
                error2_v(i) = error2_v(i) + f2(j, k)*0.1*0.1;
            end
        end
    end
end

figure
hold all
plot(mi, error2_v, 'r');
plot (mi, err2*ones(length(mi)), 'b')
xlabel('mi');
ylabel('Eps2');
legend('\epsilon2','\epsilon0');
grid on
hold off

mi = 2.245;
T = -log(mi);

figure
hold all
scatter(X1(1, :), X1(2, :), 'r*')
scatter(X2(1, :), X2(2, :), 'b*')
contour(x, y, h', [T, T], 'k', 'LineWidth', 1)
legend('Klasa 1', 'Klasa 2', 'NP klasifikator', 'Location', 'southeast')
title('Neyman-Pearson-ov klasifikator')
grid on
hold off

%% Wald sequential classifier

h1 = zeros(1, N);
for i=1:N
    X = X1(:, i);
    f11 = const1*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
    f12 = const2*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
    f21 = const3*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
    f22 = const4*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
    f1_i = P11*f11 + P12*f12;
    f2_i = P21*f21 + P22*f22;
    h1(i) = -log(f1_i/f2_i);
end

h2 = zeros(1, N);
for i=1:N
    X = X2(:, i);
    f11 = const1*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
    f12 = const2*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
    f21 = const3*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
    f22 = const4*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
    f1_i = P11*f11 + P12*f12;
    f2_i = P21*f21 + P22*f22;
    h2(i) = -log(f1_i/f2_i);
end

% definisemo vrednosti pragova za zaustavljanje algoritma
E1 = 0.000000000001;
E2 = 0.000000000001;
A = (1-E1)/E2; a = -log(A);
B = E1/(1-E2); b = -log(B);

figure(10)
hold on
plot(1:8, zeros(8), 'k')
plot(1:8, a*ones(8), 'k--')
plot(1:8, b*ones(8), 'k--')
title('Wald-ov sekvencijalni klasifikator')
xlabel('n')
ylabel('Sm')

% algoritam
for j = 1:100
    
    % racunamo Sm na osnovu slucajno izabranog oblika
    n1 = round(rand*N);
    n2 = round(rand*N);
    Sm1(1) = h1(n1);
    Sm2(1) = h2(n2);
    i = 2;
    % ukoliko nismo presli pragove, uzimamo sledeci odbirak
    while ((i <= N) && (Sm1(i-1) >= a) && (Sm1(i-1) <= b))
        Sm1(i) = Sm1(i-1) + h1(round(rand*N));
        i = i + 1;
    end

    figure(10)
    hold on
    plot(1:i-1, Sm1, 'r')
    
    i = 2;
    % ukoliko nismo presli pragove, uzimamo sledeci odbirak
    while ((i <= N) && (Sm2(i-1) >= a) && (Sm2(i-1) <= b))
        Sm2(i) = Sm2(i-1) + h2(round(rand*N));
        i = i + 1;
    end
    
    figure(10)
    hold on
    plot(1:i-1, Sm2, 'b')
    
    clear Sm1
    clear Sm2   
end

%% Dependence of m from values E1 and E2

E1 = [10^-12, 10^-10, 10^-8, 10^-6, 10^-4, 10^-3, 10^-2, 10^-1, 0.16, 0.198, 0.4, 0.98];
E2 = [10^-12, 10^-10, 10^-8, 10^-6, 10^-4, 10^-3, 10^-2, 10^-1, 0.16, 0.198, 0.4, 0.98];  
x_axis= [10^-12, 10^-10, 10^-8, 10^-6, 10^-4, 10^-3, 10^-2, 10^-1, 0.16, 0.198, 0.4, 0.98];

h1_mean = mean(h1);
m1 = (-log((1-E1)/E2(1)).*(1-E1) - log(E1/(1-E2(1))).*E1)/h1_mean;   
m2 = (-log((1-E1(1))./E2)*(1-E1(1)) - log(E1(1)./(1-E2))*E1(1))/h1_mean;

figure
hold on
semilogx(x_axis, m1)
semilogx(x_axis, m2)
ylabel('m_1');
xlabel('\xi');
title('Prosecan potreban broj odbiraka prve klase')    
legend('\xi_2 = const', '\xi_1 = const')
