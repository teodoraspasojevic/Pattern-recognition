function [N_errors] = errors_above(a, b, X)

N_errors = 0;
for i =1:length(X)
   y_real = X(2,i);
   y_est = a*X(1,i) + b;
   if y_real>y_est
       N_errors = N_errors+1;
   end
end

end

