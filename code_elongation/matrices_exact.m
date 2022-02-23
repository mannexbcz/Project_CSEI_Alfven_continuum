function [M,N] = matrices_exact(r,epsilon,k,kprime,delta,deltaprime,q,qbar,R0,B0,n,size,band,npoints)
% This function build the matrices M and N in the case of a elongated 
% and shifted equilibriuum.
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
% Output: 
%   M : matrix M, of size (size,size)
%   N : matrix N, of size (size,size)

% Useful function for the implementation of the equilibrium coefficients
    GradPsi2 = @(r,theta,k,kprime,deltaprime,qbar) B0^2.*r.^2.*k.^2.*((k.^2.*cos(theta).^2+sin(theta).^2))./(qbar.^2.*(k.*deltaprime.*cos(theta)+kprime.*r.*sin(theta).^2+k).^2); %B0^2.*r.^2.*k.^2.*
    R = @(r,theta,delta) R0+delta+r.*cos(theta);
    B2 = @(r,theta,k,kprime,delta,deltaprime,qbar) (1./R(r,theta,delta).^2).*(R0^2*B0^2+GradPsi2(r,theta,k,kprime,deltaprime,qbar));
    D = @(r,theta,k,kprime,deltaprime) kprime.*r.*sin(theta).^2+k.*deltaprime.*cos(theta)+k;
    BdotGradPhi =  @(r,theta,delta) B0*R0./(R(r,theta,delta).^2);
    BdotGradTheta = @(r,theta,k,kprime,deltaprime,delta,qbar) B0.*k./(qbar.*R(r,theta,delta).*D(r,theta,k,kprime,deltaprime));
    dThetastardTheta = @(r,theta,k,kprime,deltaprime,delta,epsilon,qbar) BdotGradPhi(r,theta,delta)./BdotGradTheta(r,theta,k,kprime,deltaprime,delta,qbar);
    F = R0*B0;
% Equilibrium coefficients
    Mtheta = @(r,theta,k,kprime,delta,deltaprime,epsilon,qbar) (B0^2/R0^2).*((R(r,theta,delta).^2)./B2(r,theta,k,kprime,delta,deltaprime,qbar)).*GradPsi2(r,theta,k,kprime,deltaprime,qbar) ;
    Ntheta = @(r,theta,k,kprime,delta,deltaprime,epsilon,qbar) GradPsi2(r,theta,k,kprime,deltaprime,qbar)./(B2(r,theta,k,kprime,delta,deltaprime,qbar).*(R(r,theta,delta).^2)) ;

    m_min=round((-n*q-size/2+1/2));
    m_max=m_min+size+1;

    m_min_coeff=0;
    m_max_coeff=(band+1)/2;
    
% Computation of theta* at middle of intervals (for midpoint quadrature)
    th=linspace(0,2*pi,npoints);
    thstar=theta_star(r,th,kprime,k,delta,deltaprime,q,qbar,B0,R0,npoints);
    thstarmid=0.5*(thstar(1:npoints-1)+thstar(2:npoints));
    
% Fourier coefficients of the equilibrium coefficients
    [coeffs_M] = get_fourier_coeff_theta_prime_exact(Mtheta,r,k,kprime,delta,deltaprime,epsilon,q,qbar,thstarmid,dThetastardTheta,m_min_coeff,m_max_coeff,npoints,R0,B0);
    [coeffs_N] = get_fourier_coeff_theta_prime_exact(Ntheta,r,k,kprime,delta,deltaprime,epsilon,q,qbar,thstarmid,dThetastardTheta,m_min_coeff,m_max_coeff,npoints,R0,B0);
       
%Construction of the matrices
    [M,N]=build_MN_bands_Fourier(coeffs_M, coeffs_N, m_min,m_max, n, q, F, size,band);

return


