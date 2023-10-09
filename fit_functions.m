function [f_total] = fit_functions(d_wall,coeffs,u_poisson)

f1 = u_poisson.^coeffs(1);
f2 = exp(-sqrt(coeffs(2))*d_wall).*(coeffs(4)*exp(sqrt(coeffs(2)-coeffs(3))* ...
    d_wall)+coeffs(5)*exp(-sqrt(coeffs(2)-coeffs(3))*d_wall));
f_total = real(f1.*f2);

end