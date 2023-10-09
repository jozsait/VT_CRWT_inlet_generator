function e = cost_function(Udata,u_poisson,coeffs,d_wall)

[fitted_data] = fit_functions(d_wall,coeffs,u_poisson);

% compute error for minimize
e = sum((fitted_data - Udata).^2);

end