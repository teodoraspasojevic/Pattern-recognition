function [V, v0] = wanted_output(X1, X2, gain1, gain2)

N = length(X1(1, :));
U = [-1*ones(1,N), ones(1,N); -1*X1, X2];
G = [ones(N,1)*gain1; ones(N, 1)*gain2];

W = inv(U*U')*U*G;

v0 = W(1);
V = [W(2); W(3)];

end

