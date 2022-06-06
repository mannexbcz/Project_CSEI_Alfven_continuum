function t = theta_star_triang(r,th,kprime,k,delta,deltaprime,d,dprime,q,qbar,B0,R0,npoints)
% This function computes the value of theta star with respect to theta,
% using equation (44)
% Inputs:
%   r : radius, double
%   th : theta, double
%   k : coefficient for the elongation, double
%   delta, deltaprime : coefficients for the shift, double
%   d,dprime : coefficients for the triangularity, double
%   q, qbar : safety factors, handle functions
%   R0 : major radius, double
%   B0 : magnetic field, double
%   npoints : number of mesh points for the numerical integration, int
% Outputs:
%   t : theta star with respect to theta, vector of size(1,lenght(th))

    Zeta = @(theta,d) theta + asin(d).*sin(theta);
    sd = @(r,d,dprime) (dprime.*r)./(sqrt(1-d.^2)) ;
    R = @(r,theta,delta,d) R0+delta+r.*cos(Zeta(theta,d));
    D = @(r,theta,k,kprime,deltaprime,d,dprime) k*cos(theta).*(deltaprime+cos(Zeta(theta,d))-sd(r,d,dprime).*sin(theta).*sin(Zeta(theta,d)))+(kprime.*r+k).*sin(theta).*(1+asin(d).*cos(theta)).*sin(Zeta(theta,d)) ;
    BdotGradPhi =  @(r,theta,delta,d) B0*R0./(R(r,theta,delta,d).^2);
    BdotGradTheta = @(r,theta,k,kprime,deltaprime,delta,d,dprime,qbar) B0.*k./(qbar.*R(r,theta,delta,d).*D(r,theta,k,kprime,deltaprime,d,dprime));

    integrand = @(theta) BdotGradPhi(r,theta,delta,d)./BdotGradTheta(r,theta,k,kprime,deltaprime,delta,d,dprime,qbar);
    
    t=zeros(1,npoints);
    sum=0;
    for i=2:npoints
        sum = sum + (1/q)*midpoint_composite_quadrature(integrand,th(i-1),th(i),1) ;
        t(i) = sum ; 
    end
 
return