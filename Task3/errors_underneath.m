function [N_errors] = errors_underneath(a, b, X)

N_errors = 0;
for i =1:length(X)
   y = X(2,i);
   y0 = a*X(1,i) + b;
   if y<y0
       N_errors = N_errors+1;
   end
end

end

