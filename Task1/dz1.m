%% Task 1 - Latter recognition

clear;
close all;
clc;

%% Loading and displayment of latters
data = load('PO_slova.mat');

N = 100;

A = data.a(:,1:N);
B = data.b(:,1:N);
C = data.c(:,1:N);
E = data.e(:,1:N);
M = data.m(:,1:N);
P = data.p(:,1:N);
S = data.s(:,1:N);
V = data.v(:,1:N);
Y = data.y(:,1:N);
Z = data.z(:,1:N);

A1 = A{1};
B1 = B{1};
C1 = C{1};
E1 = E{1};  
M1 = M{1};
P1 = P{1};
S1 = S{1};
V1 = V{1};
Y1 = Y{1};
Z1 = Z{1};

[A1x,A1y,A1p,A1xp,A1yp] = data_extraction(A1);
[B1x,B1y,B1p,B1xp,B1yp] = data_extraction(B1);
[C1x,C1y,C1p,C1xp,C1yp] = data_extraction(C1);
[E1x,E1y,E1p,E1xp,E1yp] = data_extraction(E1);
[M1x,M1y,M1p,M1xp,M1yp] = data_extraction(M1);
[P1x,P1y,P1p,P1xp,P1yp] = data_extraction(P1);
[S1x,S1y,S1p,S1xp,S1yp] = data_extraction(S1);
[V1x,V1y,V1p,V1xp,V1yp] = data_extraction(V1);
[Y1x,Y1y,Y1p,Y1xp,Y1yp] = data_extraction(Y1);
[Z1x,Z1y,Z1p,Z1xp,Z1yp] = data_extraction(Z1);

figure                 
subplot(3,2,[1 3 5]);
scatter(A1xp,A1yp)
title('Pozicija po x i y osi - A');
subplot(3,2,2)
plot(A1x);
title('Brzina po x osi - A');
subplot(3,2,4)
plot(A1y);
title('Brzina po y osi - A');
subplot(3,2,6)
plot(A1p);
title('Pritisak - A');

figure
subplot(3,2,[1 3 5]);
scatter(B1xp,B1yp)
title('Pozicija po x i y osi - B');
subplot(3,2,2)
plot(B1x);
title('Brzina po x osi - B');
subplot(3,2,4)
plot(B1y);
title('Brzina po y osi - B');
subplot(3,2,6)
plot(B1p);
title('Pritisak - B');

figure
subplot(3,2,[1 3 5]);
scatter(C1xp,C1yp)
title('Pozicija po x i y osi - C');
subplot(3,2,2)
plot(C1x);
title('Brzina po x osi - C');
subplot(3,2,4)
plot(C1y);
title('Brzina po y osi - C');
subplot(3,2,6)
plot(C1p);
title('Pritisak - C');

figure
subplot(3,2,[1 3 5]);
scatter(E1xp,E1yp)
title('Pozicija po x i y osi - E');
subplot(3,2,2)
plot(E1x);
title('Brzina po x osi - E');
subplot(3,2,4)
plot(E1y);
title('Brzina po y osi - E');
subplot(3,2,6)
plot(E1p);
title('Pritisak - E');

figure
subplot(3,2,[1 3 5]);
scatter(M1xp,M1yp)
title('Pozicija po x i y osi - M');
subplot(3,2,2)
plot(M1x);
title('Brzina po x osi - M');
subplot(3,2,4)
plot(M1y);
title('Brzina po y osi - M');
subplot(3,2,6)
plot(M1p);
title('Pritisak - M');

figure
subplot(3,2,[1 3 5]);
scatter(P1xp,P1yp)
title('Pozicija po x i y osi - P');
subplot(3,2,2)
title('Brzina po x osi - P');
plot(P1x);
subplot(3,2,4)
plot(P1y);
title('Brzina po y osi - P');
subplot(3,2,6)
plot(P1p);
title('Pritisak - P');

figure
subplot(3,2,[1 3 5]);
scatter(S1xp,S1yp)
title('Pozicija po x i y osi - S');
subplot(3,2,2)
plot(S1x);
title('Brzina po x osi - S');
subplot(3,2,4)
plot(S1y);
title('Brzina po y osi - S');
subplot(3,2,6)
plot(S1p);
title('Pritisak - S');

figure
subplot(3,2,[1 3 5]);
scatter(V1xp,V1yp)
title('Pozicija po x i y osi - V');
subplot(3,2,2)
plot(V1x);
title('Brzina po x osi - V');
subplot(3,2,4)
plot(V1y);
title('Brzina po y osi - V');
subplot(3,2,6)
plot(V1p);
title('Pritisak - V');

figure
subplot(3,2,[1 3 5]);
scatter(Y1xp,Y1yp)
title('Pozicija po x i y osi - Y');
subplot(3,2,2)
plot(Y1x);
title('Brzina po x osi - Y');
subplot(3,2,4)
plot(Y1y);
title('Brzina po y osi - Y');
subplot(3,2,6)
plot(Y1p);
title('Pritisak - Y');

figure
subplot(3,2,[1 3 5]);
scatter(Z1xp,Z1yp)
title('Pozicija po x i y osi - Z');
subplot(3,2,2)
plot(Z1x);
title('Brzina po x osi - Z');
subplot(3,2,4)
plot(Z1y);
title('Brzina po y osi - Z');
subplot(3,2,6)
plot(Z1p);
title('Pritisak - Z');

%% Faeture extraction

N_features = 15;
FA = zeros(N_features, N);
FB = zeros(N_features, N);
FC = zeros(N_features, N);
FE = zeros(N_features, N);
FM = zeros(N_features, N);
FP = zeros(N_features, N);
FS = zeros(N_features, N);
FV = zeros(N_features, N);
FY = zeros(N_features, N);
FZ = zeros(N_features, N);

for i=1:N
   FA(:, i) = feature_extraction(A{i});
   FB(:, i) = feature_extraction(B{i});
   FC(:, i) = feature_extraction(C{i});
   FE(:, i) = feature_extraction(E{i});
   FM(:, i) = feature_extraction(M{i});
   FP(:, i) = feature_extraction(P{i});
   FS(:, i) = feature_extraction(S{i});
   FV(:, i) = feature_extraction(V{i});
   FY(:, i) = feature_extraction(Y{i});
   FZ(:, i) = feature_extraction(Z{i}); 
end

%% Dimensionality reduction uding LTA method

MA = mean(FA, 2); SA = cov(FA');
MB = mean(FB, 2); SB = cov(FB');
MC = mean(FC, 2); SC = cov(FC');
ME = mean(FE, 2); SE = cov(FE');
MM = mean(FM, 2); SM = cov(FM');
MP = mean(FP, 2); SP = cov(FP');
MS = mean(FS, 2); SS = cov(FS');
MV = mean(FV, 2); SV = cov(FV');
MY = mean(FY, 2); SY = cov(FY');
MZ = mean(FZ, 2); SZ = cov(FZ');

Pi = 1/10;

Sw = (SA + SB + SC + SE + SM + SP + SS + SV + SY + SZ)/Pi;

M0 = (MA + MB + MC + ME + MM + MP + MS + MV + MY + MZ)/Pi;
Sb = ((MA-M0)*(MA-M0)' + (MB-M0)*(MB-M0)' + (MC-M0)*(MC-M0)' +...
      (ME-M0)*(ME-M0)' + (MM-M0)*(MM-M0)' + (MP-M0)*(MP-M0)' + ...
      (MS-M0)*(MS-M0)' + (MY-M0)*(MY-M0)' + (MV-M0)*(MV-M0)' + ...
      (MZ-M0)*(MZ-M0)')/Pi;

S = Sw^(-1)*Sb;

[F, L] = eig(S);
transformator(:,1) = F(:,8);
transformator(:,2) = F(:,9);
transformator(:,3) = F(:,2);

FA_red = transformator'*FA;
FB_red = transformator'*FB;
FC_red = transformator'*FC;
FE_red = transformator'*FE;
FM_red = transformator'*FM;
FP_red = transformator'*FP;
FS_red = transformator'*FS;
FV_red = transformator'*FV;
FY_red = transformator'*FY;
FZ_red = transformator'*FZ;

figure
hold all
scatter3(FA_red(1,:),FA_red(2,:),FA_red(3,:),'r*');
scatter3(FB_red(1,:),FB_red(2,:),FB_red(3,:),'g*');
scatter3(FC_red(1,:),FC_red(2,:),FC_red(3,:),'b*');
scatter3(FE_red(1,:),FE_red(2,:),FE_red(3,:),'c*');
scatter3(FM_red(1,:),FM_red(2,:),FM_red(3,:),'m*');
scatter3(FP_red(1,:),FP_red(2,:),FP_red(3,:),'y*');
scatter3(FS_red(1,:),FS_red(2,:),FS_red(3,:),'k*');
scatter3(FV_red(1,:),FV_red(2,:),FV_red(3,:),'ro');
scatter3(FY_red(1,:),FY_red(2,:),FY_red(3,:),'go');
scatter3(FZ_red(1,:),FZ_red(2,:),FZ_red(3,:),'bo');
legend('A','B','C','E','M','P','S','V','Y','Z')
title('Obelezja slova')
grid on
hold off
 
%% Classificator design

MA = mean(FA_red, 2); SA = cov(FA_red');
MB = mean(FB_red, 2); SB = cov(FB_red');
MC = mean(FC_red, 2); SC = cov(FC_red');
ME = mean(FE_red, 2); SE = cov(FE_red');
MM = mean(FM_red, 2); SM = cov(FM_red');
MP = mean(FP_red, 2); SP = cov(FP_red');
MS = mean(FS_red, 2); SS = cov(FS_red');
MV = mean(FV_red, 2); SV = cov(FV_red');
MY = mean(FY_red, 2); SY = cov(FY_red');
MZ = mean(FZ_red, 2); SZ = cov(FZ_red');

classes = ['A','B','C','E','M','P','S','V','Y','Z'];

conf_matrix = zeros(10,10);

for k = 1:length(classes)
   if classes(k) == 'A' 
       F = FA_red;
   elseif classes(k) == 'B'
       F = FB_red;
   elseif classes(k) == 'C'
       F = FC_red;
   elseif classes(k) == 'E'
       F = FE_red;
   elseif classes(k) == 'M'
       F = FM_red;
   elseif classes(k) == 'P'
       F = FP_red;
   elseif classes(k) == 'S'
       F = FS_red;
   elseif classes(k) == 'V'
       F = FV_red;
   elseif classes(k) == 'Y'
       F = FY_red;
   elseif classes(k) == 'Z'
       F = FZ_red;
   end
    
   for i =1:length(F(1,:))
        X = F(:,i);
        fa = 1/((2*pi)^(3/2)*det(SA)^0.5)*exp(-1/2*(X-MA)'*inv(SA)*(X-MA));
        fb = 1/((2*pi)^(3/2)*det(SB)^0.5)*exp(-1/2*(X-MB)'*inv(SA)*(X-MB));
        fc = 1/((2*pi)^(3/2)*det(SC)^0.5)*exp(-1/2*(X-MC)'*inv(SC)*(X-MC));
        fe = 1/((2*pi)^(3/2)*det(SE)^0.5)*exp(-1/2*(X-ME)'*inv(SE)*(X-ME));
        fm = 1/((2*pi)^(3/2)*det(SM)^0.5)*exp(-1/2*(X-MM)'*inv(SM)*(X-MM));
        fp = 1/((2*pi)^(3/2)*det(SP)^0.5)*exp(-1/2*(X-MP)'*inv(SP)*(X-MP));
        fs = 1/((2*pi)^(3/2)*det(SS)^0.5)*exp(-1/2*(X-MS)'*inv(SS)*(X-MS));
        fv = 1/((2*pi)^(3/2)*det(SV)^0.5)*exp(-1/2*(X-MV)'*inv(SV)*(X-MV));
        fy = 1/((2*pi)^(3/2)*det(SY)^0.5)*exp(-1/2*(X-MY)'*inv(SY)*(X-MY));
        fz = 1/((2*pi)^(3/2)*det(SZ)^0.5)*exp(-1/2*(X-MZ)'*inv(SZ)*(X-MZ));
        
        ml = max([fa,fb,fc,fe,fm,fp,fs,fv,fy,fz]);
        
        if ml == fa
            conf_matrix(1,k)= conf_matrix(1,k)+1;
        elseif ml == fb
            conf_matrix(2,k)= conf_matrix(2,k)+1;
        elseif ml == fc
            conf_matrix(3,k)= conf_matrix(3,k)+1; 
        elseif ml == fe
            conf_matrix(4,k)= conf_matrix(4,k)+1;
        elseif ml == fm
            conf_matrix(5,k)= conf_matrix(5,k)+1;
        elseif ml == fp
            conf_matrix(6,k)= conf_matrix(6,k)+1;
        elseif ml == fs
            conf_matrix(7,k)= conf_matrix(7,k)+1;
        elseif ml == fv
            conf_matrix(8,k)= conf_matrix(8,k)+1;
        elseif ml == fy
            conf_matrix(9,k)= conf_matrix(9,k)+1;
        elseif ml == fz
            conf_matrix(10,k)= conf_matrix(10,k)+1;
        end
    end
end

disp('Konfuziona matrica:');
disp(conf_matrix)

error = (sum(sum(conf_matrix))-trace(conf_matrix))/sum(sum(conf_matrix));
text = ['Ukupna greska: ', num2str(error), '%'];
disp(text);

%% Correctly and wrongly classified latter M

max_likelihood_e = 0;
max_likelihood_m = 0;

for i =1:length(FM_red(1,:))
        X = FM_red(:,i);
        fa = 1/((2*pi)^(3/2)*det(SA)^0.5)*exp(-1/2*(X-MA)'*inv(SA)*(X-MA));
        fb = 1/((2*pi)^(3/2)*det(SB)^0.5)*exp(-1/2*(X-MB)'*inv(SA)*(X-MB));
        fc = 1/((2*pi)^(3/2)*det(SC)^0.5)*exp(-1/2*(X-MC)'*inv(SC)*(X-MC));
        fe = 1/((2*pi)^(3/2)*det(SE)^0.5)*exp(-1/2*(X-ME)'*inv(SE)*(X-ME));
        fm = 1/((2*pi)^(3/2)*det(SM)^0.5)*exp(-1/2*(X-MM)'*inv(SM)*(X-MM));
        fp = 1/((2*pi)^(3/2)*det(SP)^0.5)*exp(-1/2*(X-MP)'*inv(SP)*(X-MP));
        fs = 1/((2*pi)^(3/2)*det(SS)^0.5)*exp(-1/2*(X-MS)'*inv(SS)*(X-MS));
        fv = 1/((2*pi)^(3/2)*det(SV)^0.5)*exp(-1/2*(X-MV)'*inv(SV)*(X-MV));
        fy = 1/((2*pi)^(3/2)*det(SY)^0.5)*exp(-1/2*(X-MY)'*inv(SY)*(X-MY));
        fz = 1/((2*pi)^(3/2)*det(SZ)^0.5)*exp(-1/2*(X-MZ)'*inv(SZ)*(X-MZ));
        
        ml = max([fa,fb,fc,fe,fm,fp,fs,fv,fy,fz]);
        
        if ml == fe
            if fe > max_likelihood_e
                max_likelihood_e = fe;
                index_as_E = i;
            end
        elseif ml == fm
            if fm > max_likelihood_m
                max_likelihood_m = fm;
                index_as_M = i;
            end
        end
end

[Mwx,Mwy,Mwp,Mwxp,Mwyp] = data_extraction(M{index_as_E});
[Mcx,Mcy,Mcp,Mcxp,Mcyp] = data_extraction(M{index_as_M});

figure
scatter(Mwxp,Mwyp)
title('Slovo M, pogresno klasifikovano kao E');

figure
scatter(Mcxp,Mcyp)
title('Slovo M, tacno klasifikovano kao M');

%% Classification of two latters using parametric method

% latters E and S

gainE = 3;
gainS = 1;

[V, v0] = wanted_output(FE_red, FS_red, gainE, gainS);

x = linspace(-6,6,N);
y = linspace(-6,8,N);
classificator = V'*[x; y] + v0;

figure
hold all
scatter3(FE_red(1,:),FE_red(2,:),FE_red(3,:),'r*');
scatter3(FS_red(1,:),FS_red(2,:),FS_red(3,:),'b*');
scatter3(x, y, classificator, 'k.')
legend('Obelezja slova E', 'Obelezja slova S', 'Klasifikator')
title('Klasifikator dobijen metodom zeljenog izlaza')
grid on
hold off
