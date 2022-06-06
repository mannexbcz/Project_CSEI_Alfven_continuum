function [coeffs]=get_fourier_coeff_chease(fun,thetastar,mmin,mmax)
% This function returns the Fourier coefficients betweem nmin and nmax of
% fun
% Inputs:
%       fun : vector of doubles
%       thetastar : vector of doubles
%       dThetastardTheta : partial derivative of theta star with respect to
%                           theta, handle function
%       mmin,mmax : limiting mode number to be computed
%       R0 : major radius, double
%       B0 : magnetic fiels, double
%       npoints : number of mesh points for the numerical integration, int
% Output:
%   coeffs : fourier coefficients, vector of doubles
   
    coeffs=[];
    thetastarlong = [thetastar,2*pi];
    for m = mmin:mmax
        fun_theta = exp(-1i.*m.*thetastar).*fun;
        % Interpolation of the coefficient in thetastar = 2*pi to ensure an
        % integration between 0 and 2*pi
        fit = polyfit(thetastar(end-1:end), fun_theta(end-1:end),1);
        fun_end = fun_theta(1);
        fun_theta=[fun_theta,fun_end];
        % Integration via trapezoidal rule
        c = (1/(2*pi)).*trapz(thetastarlong,fun_theta);
        coeffs=[coeffs,c];
        
    end
end