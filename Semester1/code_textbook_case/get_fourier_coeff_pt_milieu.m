function [coeffs]=get_fourier_coeff_pt_milieu(fun,epsilon,epsilon2,epsilon3,nmin,nmax,npoints)
% This function returns the Fourier coefficients betweem nmin and nmax of
% the function fun
% Inputs:
%       fun : handle function
%       epsilon, epsilon2, epsilon3 : coefficients used in the function fun, doubles
%       nmin,nmax : limiting indices of the Fourier coefficients computed, int
%       npoints : number of mesh points for the numerical integration, int
% Output:
%       coeffs : Fourier coefficients

    coeffs=[];
    for m = nmin:nmax
        fun_theta = @(theta) exp(-1i.*m.*theta).*fun(epsilon,epsilon2,epsilon3,theta); 
        c=(1/2*pi)*midpoint_composite_quadrature(fun_theta, 0, 2*pi, npoints);       
        coeffs=[coeffs,c];
    end

end