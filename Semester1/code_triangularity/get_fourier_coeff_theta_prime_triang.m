function [coeffs]=get_fourier_coeff_theta_prime_triang(fun,r,k,kprime,delta,deltaprime,d,dprime,epsilon,q,qbar,thetastar,dThetastardTheta,mmin,mmax,npoints,R0,B0)
% This function returns the Fourier coefficients betweem nmin and nmax of
% the function fun
% Inputs:
%       fun : handle function
%       r : radius, double
%       k, kprime : coefficients for the elongation, double
%       delta, deltaprime : coefficients for the shift, double
%       d,dprime : coefficients for the triangularity
%       q, qbar : safety factors, doubles
%       thetastar : vector
%       dThetastardTheta : partial derivative of theta star with respect to
%                           theta, handle function
%       mmin,mmax : limiting mode number to be computed
%       R0 : major radius, double
%       B0 : magnetic fiels, double
%       npoints : number of mesh points for the numerical integration, int
% Output:
%   coeffs : fourier coefficients, vector of doubles
   
    coeffs=[];
    for m = mmin:mmax
        fun_theta = @(theta) exp(-1i.*m.*thetastar).*fun(r,theta,k,kprime,delta,deltaprime,d,dprime,epsilon,qbar).*dThetastardTheta(r,theta,k,kprime,deltaprime,delta,d,dprime,epsilon,qbar);
        c=(1/2*pi)*midpoint_composite_quadrature(fun_theta, 0, 2*pi, npoints-1);
        coeffs=[coeffs,c];
    end

end