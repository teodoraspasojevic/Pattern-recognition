function [Y1, Y2, l] = sq_decomp(Y1, Y2)

N1 = length(Y1(1, :));
N2 = length(Y2(1, :));
N = N1 + N2;

M1 = mean(Y1, 2);
M2 = mean(Y2, 2);

S1 = cov(Y1');
S2 = cov(Y2');

P1 = N1/N;
P2 = N2/N;

lmax = 100;
l = 0;
reclasterization = 1;
while((l<lmax) && reclasterization)
    
    l = l + 1;
    
    Y1_pom = [];
    Y2_pom = [];
    reclasterization = 0;
    
    if N1~=0
        for i=1:N1
            d1 = 0.5*(Y1(:,i)-M1)'*S1^(-1)*(Y1(:,i)-M1) + 0.5*log(det(S1))-0.5*log(P1);
            d2 = 0.5*(Y1(:,i)-M2)'*S2^(-1)*(Y1(:,i)-M2) + 0.5*log(det(S2))-0.5*log(P2);
            d = [d1 d2];
            d_min = min(d);
            if d_min == d1
                Y1_pom = [Y1_pom Y1(:, i)];
            else
                Y2_pom = [Y2_pom Y1(:, i)];
            end
            if d_min ~= d1
                reclasterization = 1;
            end
        end
    end
    if N2~=0
        for i=1:N2
            d1 = 0.5*(Y2(:,i)-M1)'*S1^(-1)*(Y2(:,i)-M1) + 0.5*log(det(S1))-0.5*log(P1);
            d2 = 0.5*(Y2(:,i)-M2)'*S2^(-1)*(Y2(:,i)-M2) + 0.5*log(det(S2))-0.5*log(P2);
            d = [d1 d2];
            d_min = min(d);
            if d_min == d1
                Y1_pom = [Y1_pom Y2(:, i)];
            else
                Y2_pom = [Y2_pom Y2(:, i)];
            end
            if d_min ~= d2
                reclasterization = 1;
            end
        end
    end
    
    clear Y1 Y2
    Y1 = Y1_pom;
    Y2 = Y2_pom;
    
    if isempty(Y1)
        N1 = 0;
        M1 = Inf;
    else
        N1 = length(Y1(1, :));
        M1 = mean(Y1, 2);
    end
    if isempty(Y2)
        N2 = 0;
        M1 = Inf;
    else
        N2 = length(Y2(1, :));
        M2 = mean(Y2, 2);
    end
    N = N1 + N2;

    S1 = cov(Y1');
    S2 = cov(Y2');

    P1 = N1/N;
    P2 = N2/N;
    
end

text = ['Potreban broj iteracija je ', num2str(l), '.' ];
disp(text);


figure
hold all
if ~isempty(Y1)
    scatter(Y1(1, :), Y1(2, :), 'r*')
end
if ~isempty(Y2)
    scatter(Y2(1, :), Y2(2, :), 'b*')
end
title('Odbirci klasa posle kvadratne dekompozicije')
grid on
hold off

end

