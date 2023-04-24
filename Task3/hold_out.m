function [V_opt, v0_opt] = hold_out(X1_train, X2_train, X1_test, X2_test)

M1_train = mean(X1_train, 2);
M2_train = mean(X2_train, 2);

S1_train = cov(X1_train');
S2_train = cov(X2_train');

% na train skupu estimiramo parametre linearnog klasifikatora - v0, V
% na test skupu brojimo pogresno klasifikovane oblike

s = 0:0.001:1;

% u 1. vrsi cuvamo broj gresaka u datoj iteraciji
% u 2. vrsti cuvamo v0 za datu iteraciju
% u 3. vrsti cuvamo s u datoj iteraciji
parameters = zeros(3, length(s));

for i=1:length(s)
   
    % 1. nalazenje vektora V
    V = (s(i)*S1_train + (1-s(i))*S2_train)^(-1)*(M2_train - M1_train);
    
    % 2. projektovanje vektora X(test) na vektor V, za brojanje gresaka
    Y1err = V'*X1_test;
    Y2err = V'*X2_test;
    
    % 3. projektovanje vektora X(train) na vektor V, za nalazenje praga -v0
    Y1 = V'*X1_train;
    Y2 = V'*X2_train;
    
    Y = sort([Y1 Y2]);
    
    % 4. brojenje gresaka (odbirci prokejtovani na V kao Y1, a v0 ih
    % klasifikuje kao Y2 i odbirci projektovani na V kao Y2, a v0 ih
    % klasifikuje kao Y1
    v0 = zeros(1,length(Y)-1);
    N_error = zeros(1,length(Y)-1);
    
    for j=1:length(v0)
        
        v0(j) = -(Y(j) + Y(j+1))/2;
        
        for k=1:length(Y1err)
            if Y1err(k) > -v0(j)
                N_error(j) = N_error(j) + 1;
            end
        end
        for k=1:length(Y2err)
            if Y2err(k) < -v0(j)
                N_error(j) = N_error(j) + 1;
            end
        end
    end
    
    % 5. nalazenje optimalnog v0 u zavisnosti od N_error
    [N_error_min, index] = min(N_error);
    v0_opt = v0(index);
    
    % 6. cuvanje vrednosti N_error, v0 i s po iteracijama
    parameters(1, i) = N_error_min;
    parameters(2, i) = v0_opt;
    parameters(3, i) = s(i);
    
end

[N_error_min, index] = min(parameters(1, :));
v0_opt = parameters(2, index);
s_opt = parameters(3, index);
V_opt =  (s_opt*S1_train + (1-s_opt)*S2_train)^(-1)*(M2_train - M1_train);

end

