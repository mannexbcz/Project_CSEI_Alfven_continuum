function [coeffs]=get_fourier_coeff_chease(fun,thetastar,mmin,mmax)
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
        fun_theta = exp(-1i.*m.*thetastar).*fun;
        thmid = (thetastar(2:end)-thetastar(1:end-1))./2;
        intervals = zeros(1,length(fun));
        intervals(2:end-1)=thmid(2:end)-thmid(1:end-1);
        intervals(1)=thmid(1);
        intervals(end)=2*pi-thmid(end);
        c=(1/(2*pi)).*intervals.*fun_theta;
        coeffs=[coeffs,c];
        
    end
    disp(length(coeffs))
save('coeffs.mat', 'coeffs')
end