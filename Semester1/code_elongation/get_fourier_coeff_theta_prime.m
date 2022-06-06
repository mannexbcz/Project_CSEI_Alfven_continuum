function [coeffs]=get_fourier_coeff_theta_prime(fun,r,k,kprime,delta,deltaprime,epsilon,qbar,thetastar,dThetastardTheta,mmin,mmax,npoints)
% This function returns the Fourier coefficients betweem nmin and nmax of
% the function fun
% Inputs:
%       fun : handle function
%       r : radius, double
%       k, kprime : coefficients for the elongation, double
%       delta, deltaprime : coefficients for the shift, double
%       qbar : safety factor, double
%       thetastar : handle function giving thetastar wrt theta
%       dThetastardTheta : partial derivative of theta star with respect to
%                           theta, handle function
%       mmin,mmax : limiting mode number to be computed
%       npoints : number of mesh points for the numerical integration, int
    
    coeffs=[];
    for m = mmin:mmax
        fun_theta = @(theta) exp(-1i.*m.*thetastar(r,theta,k,kprime,deltaprime,epsilon)).*fun(r,theta,k,kprime,delta,deltaprime,epsilon,qbar).*dThetastardTheta(r,theta,k,kprime,deltaprime,delta,epsilon,qbar);
        c=(1/2*pi)*midpoint_composite_quadrature(fun_theta, 0, 2*pi, npoints-1);
        coeffs=[coeffs,c];
    end

end