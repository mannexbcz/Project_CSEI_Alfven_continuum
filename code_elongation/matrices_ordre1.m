function [M,N] = matrices_ordre1(r,epsilon,k,kprime,delta,deltaprime,q,qbar,R0,B0,n,size,band,npoints)
% This function build the matrices M and N in the case of a elongated 
% and shifted equilibriuum, and in the first order approximation. The
% equilibrium coefficients are given by equations (62) and (63)
% Inputs:
%   r : radius, double
%   epsilon : large aspect ratio, double
%   k, kprime : coefficients for the elongation, double
%   delta, deltaprime : coefficients for the shift, double
%   q, qbar : safety factors, doubles
%   n : toroidal mode number, int
%   R0 : major radius, double
%   B0 : magnetic fiels, double
%   size : required rank of the matrices, int
%   band : required number of bands of the matrices, int
%   npoints : number of mesh points for the numerical integration, int
% Outputs: 
%   M : matrix M, of size (size,size)
%   N : matrix N, of size (size,size)

F=R0*B0;
% Useful functions
thetastar = @(r,theta,k,kprime,deltaprime,epsilon) theta+(deltaprime-epsilon).*sin(theta)+(kprime.*r./(4.*k))*sin(2.*theta);
dThetastardTheta = @(r,theta,k,kprime,deltaprime,delta,epsilon,qbar) (kprime.*r./k).*sin(theta).^2+(deltaprime-epsilon).*cos(theta)+1;

% Equilibrium coefficients
Mtheta = @(r,theta,k,kprime,delta,deltaprime,epsilon,qbar) (B0^4*R0^4)./(qbar.^2).*(k.^2.*cos(theta).^2+sin(theta).^2).*(1-(2*kprime*R0/k).*epsilon.*sin(theta).^2-(2.*deltaprime-4.*epsilon).*cos(theta));
Ntheta = @(r,theta,k,kprime,delta,deltaprime,epsilon,qbar) (B0^2*R0^2)./(qbar.^2).*(k.^2.*cos(theta).^2+sin(theta).^2).*(1-(2*kprime*R0/k).*epsilon.*sin(theta).^2-2.*deltaprime.*cos(theta));

    m_min=round((-n*q-size/2+1/2));
    m_max=m_min+size+1;
    m_min_coeff=0;
    m_max_coeff=(band+1)/2;
      
    [coeffs_M] = get_fourier_coeff_theta_prime(Mtheta,r,k,kprime,delta,deltaprime,epsilon,qbar,thetastar,dThetastardTheta,m_min_coeff,m_max_coeff,npoints);
    [coeffs_N] = get_fourier_coeff_theta_prime(Ntheta,r,k,kprime,delta,deltaprime,epsilon,qbar,thetastar,dThetastardTheta,m_min_coeff,m_max_coeff,npoints);
       
    %Construction of the matrices
    [M,N]=build_MN_bands_Fourier(coeffs_M, coeffs_N, m_min,m_max, n, qbar, F, size,band);
return