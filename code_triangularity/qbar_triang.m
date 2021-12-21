function qbar = qbar_triang(q,r,k,kprime,delta,deltaprime,d,dprime,R0,npoints)
% This function implements the numerical computation of qbar according to
% equation (21) and for the case of an elongated and shifted equilibriuum,
% with triangularity.
% Inputs:
%       q : safety factor, handle function
%       r : radius, double
%       k, kprime : coefficients for the elongation, double
%       delta, deltaprime : coefficients for the shift, double
%       d,dprime : coefficients for the triangularity, double
%       R0 : major radius, double
%       npoints : number of mesh points for the numerical integration, int
% Outputs:
%       qbar, double

    Zeta = @(theta,d) theta + asin(d).*sin(theta);
    R = @(r,theta,delta,d) R0+delta+r.*cos(Zeta(theta,d));
    sd = @(r,d,dprime) (dprime.*r)./(sqrt(1-d.^2)) ;
    D = @(r,theta,k,kprime,deltaprime,d,dprime) k*cos(theta).*(deltaprime+cos(Zeta(theta,d))+sd(r,d,dprime).*sin(theta).*sin(Zeta(theta,d)))+(kprime.*r+k).*sin(theta).*(1+asin(d).*cos(theta)).*sin(Zeta(theta,d)) ;

    integrand_q = @(theta)  (R0.*D(r,theta,k,kprime,deltaprime,d,dprime))./(R(r,theta,delta,d).*k) ;
    
    qbar = q*2*pi/midpoint_composite_quadrature(integrand_q, 0, 2*pi, npoints);
end