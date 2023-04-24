function [Y1, Y2, Y3, Y4, l] = c_mean(Y1, Y2, Y3, Y4)

N1 = length(Y1(1, :));
N2 = length(Y2(1, :));
N3 = length(Y3(1, :));
N4 = length(Y4(1, :));

M1 = mean(Y1, 2);
M2 = mean(Y2, 2);
M3 = mean(Y3, 2);
M4 = mean(Y4, 2);

lmax = 100;
l = 0;
reclasterization = 1;
while((l<lmax) && reclasterization)
    
    l = l + 1;
    
    Y1_pom = [];
    Y2_pom = [];
    Y3_pom = [];
    Y4_pom = [];
    reclasterization = 0;
    
    if N1~=0
        for i=1:N1
            d1 = sum((Y1(:,i)-M1).^2);
            d2 = sum((Y1(:,i)-M2).^2);
            d3 = sum((Y1(:,i)-M3).^2);
            d4 = sum((Y1(:,i)-M4).^2);
            d = [d1 d2 d3 d4];
            d_min = min(d);
            if d_min == d1
                Y1_pom = [Y1_pom Y1(:, i)];
            elseif d_min == d2
                Y2_pom = [Y2_pom Y1(:, i)];
            elseif d_min == d3
                Y3_pom = [Y3_pom Y1(:, i)];
            else
                Y4_pom = [Y4_pom Y1(:, i)];
            end
            if d_min ~= d1
                reclasterization = 1;
            end
        end
    end
    if N2~=0
        for i=1:N2
            d1 = sum((Y2(:,i)-M1).^2);
            d2 = sum((Y2(:,i)-M2).^2);
            d3 = sum((Y2(:,i)-M3).^2);
            d4 = sum((Y2(:,i)-M4).^2);
            d = [d1 d2 d3 d4];
            d_min = min(d);
            if d_min == d1
                Y1_pom = [Y1_pom Y2(:, i)];
            elseif d_min == d2
                Y2_pom = [Y2_pom Y2(:, i)];
            elseif d_min == d3
                Y3_pom = [Y3_pom Y2(:, i)];
            else
                Y4_pom = [Y4_pom Y2(:, i)];
            end
            if d_min ~= d2
                reclasterization = 1;
            end
        end
    end
    if N3~=0
        for i=1:N3
            d1 = sum((Y3(:,i)-M1).^2);
            d2 = sum((Y3(:,i)-M2).^2);
            d3 = sum((Y3(:,i)-M3).^2);
            d4 = sum((Y3(:,i)-M4).^2);
            d = [d1 d2 d3 d4];
            d_min = min(d);
            if d_min == d1
                Y1_pom = [Y1_pom Y3(:, i)];
            elseif d_min == d2
                Y2_pom = [Y2_pom Y3(:, i)];
            elseif d_min == d3
                Y3_pom = [Y3_pom Y3(:, i)];
            else
                Y4_pom = [Y4_pom Y3(:, i)];
            end
            if d_min ~= d3
                reclasterization = 1;
            end
        end
    end
    if N4~=0
        for i=1:N4
            d1 = sum((Y4(:,i)-M1).^2);
            d2 = sum((Y4(:,i)-M2).^2);
            d3 = sum((Y4(:,i)-M3).^2);
            d4 = sum((Y4(:,i)-M4).^2);
            d = [d1 d2 d3 d4];
            d_min = min(d);
            if d_min == d1
                Y1_pom = [Y1_pom Y4(:, i)];
            elseif d_min == d2
                Y2_pom = [Y2_pom Y4(:, i)];
            elseif d_min == d3
                Y3_pom = [Y3_pom Y4(:, i)];
            else
                Y4_pom = [Y4_pom Y4(:, i)];
            end
            if d_min ~= d4
                reclasterization = 1;
            end 
        end
    end
    
    clear Y1 Y2 Y3 Y4
    Y1 = Y1_pom;
    Y2 = Y2_pom;
    Y3 = Y3_pom;
    Y4 = Y4_pom;
    
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
    if isempty(Y3)
        N3 = 0;
        M1 = Inf;
    else
        N3 = length(Y3(1, :));
        M3 = mean(Y3, 2);
    end
    if isempty(Y4)
        N4 = 0;
        M1 = Inf;
    else
        N4 = length(Y4(1, :));
        M4 = mean(Y4, 2);
    end
    
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
if ~isempty(Y3)
    scatter(Y3(1, :), Y3(2, :), 'g*')
end
if ~isempty(Y4)
    scatter(Y4(1, :), Y4(2, :), 'y*')
end
title('Odbirci klasa posle c mean klasterizacije')
grid on
hold off

end

