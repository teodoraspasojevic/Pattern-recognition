function [Q, V, v0] = wanted_output_sq(X1, X2, gain1, gain2)

N = length(X1(1, :));
U = [-1*ones(1,N), ones(1,N); -1*X1, X2; -1*(X1(1,:)).^2,(X2(1,:)).^2;...
    -1*(X1(2,:)).^2, (X2(2,:)).^2; -2*X1(1,:).*X1(2,:), 2*X2(1,:).*X2(2,:)];
G = [ones(N,1)*gain1; ones(N, 1)*gain2];

W = inv(U*U')*U*G;

v0 = W(1);
V = [W(2); W(3)];
Q = [W(4); W(5); W(6)];

end

