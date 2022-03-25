function [M,N] = matrices_triang(r,epsilon,k,kprime,delta,deltaprime,d,dprime,q,qbar,R0,B0,n,size,band,npoints)
% This function build the matrices M and N in the case of a elongated 
% and shifted equilibriuum, with triangularity.
% Inputs:
%   r : radius, double
%   epsilon : large aspect ratio, double
%   k, kprime : coefficients for the elongation, double
%   delta, deltaprime : coefficients for the shift, double
%   d,dprime : coefficients for the triangularity
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

% Useful functions
Zeta = @(theta,d) theta + asin(d).*sin(theta);
sd = @(r,d,dprime) (dprime.*r)./(sqrt(1-d.^2)) ;
sk = @(r,k,kprime) (kprime.*r)./k ; 
R = @(r,theta,delta,d) R0+delta+r.*cos(Zeta(theta,d));
D = @(r,theta,k,kprime,deltaprime,d,dprime) k*cos(theta).*(deltaprime+cos(Zeta(theta,d))-sd(r,d,dprime).*sin(theta).*sin(Zeta(theta,d)))+(kprime.*r+k).*sin(theta).*(1+asin(d).*cos(theta)).*sin(Zeta(theta,d)) ;
GradPsi2 = @(r,theta,k,kprime,deltaprime,d,dprime,qbar) B0^2.*r.^2.*k.^2.*(k.^2.*cos(theta).^2+(1+asin(d).*cos(theta)).^2.*sin(Zeta(theta,d)).^2)./(qbar.^2.*D(r,theta,k,kprime,deltaprime,d,dprime).^2); 
B2 = @(r,theta,k,kprime,delta,deltaprime,d,dprime,qbar) (1./R(r,theta,delta,d).^2).*(R0^2*B0^2+GradPsi2(r,theta,k,kprime,deltaprime,d,dprime,qbar));
BdotGradPhi =  @(r,theta,delta,d) B0*R0./(R(r,theta,delta,d).^2);
BdotGradTheta = @(r,theta,k,kprime,deltaprime,delta,d,dprime,qbar) B0.*k./(qbar.*R(r,theta,delta,d).*D(r,theta,k,kprime,deltaprime,d,dprime));
dThetastardTheta = @(r,theta,k,kprime,deltaprime,delta,d,dprime,epsilon,qbar) BdotGradPhi(r,theta,delta,d)./BdotGradTheta(r,theta,k,kprime,deltaprime,delta,d,dprime,qbar);
F = R0*B0;
% Equilibrium coefficients
Mtheta = @(r,theta,k,kprime,delta,deltaprime,d,dprime,epsilon,qbar) (B0^2/R0^2).*((R(r,theta,delta,d).^2)./B2(r,theta,k,kprime,delta,deltaprime,d,dprime,qbar)).*GradPsi2(r,theta,k,kprime,deltaprime,d,dprime,qbar) ;
Ntheta = @(r,theta,k,kprime,delta,deltaprime,d,dprime,epsilon,qbar) GradPsi2(r,theta,k,kprime,deltaprime,d,dprime,qbar)./(B2(r,theta,k,kprime,delta,deltaprime,d,dprime,qbar).*(R(r,theta,delta,d).^2)) ;

    m_min=round((-n*q-size/2+1/2));
    m_max=m_min+size+1;

    m_min_coeff=0;
    m_max_coeff=(band+1)/2;
    
   % Computation of theta * at middle of intervals (for midpoint quadrature)
    th=linspace(0,2*pi,npoints);
    thstar=theta_star_triang(r,th,kprime,k,delta,deltaprime,d,dprime,q,qbar,B0,R0,npoints);
    thstarmid=0.5*(thstar(1:npoints-1)+thstar(2:npoints));
%     thmid=0.5*(th(1:end-1)+th(2:end));
%     thstarmid=theta_star_triang(r,thmid,kprime,k,delta,deltaprime,d,dprime,q,qbar,B0,R0,npoints-1);
    
%     thmid = 0.5*(th(1:end-1)+th(2:end));
%     thstarmid = theta_star_mid(r,th,kprime,k,delta,deltaprime,d,dprime,q,qbar,B0,R0,npoints-1);

   % Fourier coefficients
    [coeffs_M] = get_fourier_coeff_theta_prime_triang(Mtheta,r,k,kprime,delta,deltaprime,d,dprime,epsilon,q,qbar,thstarmid,dThetastardTheta,m_min_coeff,m_max_coeff,npoints,R0,B0);
    [coeffs_N] = get_fourier_coeff_theta_prime_triang(Ntheta,r,k,kprime,delta,deltaprime,d,dprime,epsilon,q,qbar,thstarmid,dThetastardTheta,m_min_coeff,m_max_coeff,npoints,R0,B0);

   %Construction of the matrices
    [M,N]=build_MN_bands_Fourier(coeffs_M, coeffs_N, m_min,m_max, n, q, F, size,band);
return


