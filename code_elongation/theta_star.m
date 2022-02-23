function t = theta_star(r,th,kprime,k,delta,deltaprime,q,qbar,B0,R0,npoints)
% This function computes the value of theta star with respect to theta,
% using equation (44)
% Inputs:
%   r : radius, double
%   th : theta, double
%   k, kprime : coefficients for the elongation, double
%   delta, deltaprime : coefficients for the shift, double
%   q, qbar : safety factors, handle functions
%   R0 : major radius, double
%   B0 : magnetic field, double
%   npoints : number of mesh points for the numerical integration, int
% Outputs:
%   t : theta star with respect to theta, vector of size(1,lenght(th))

    R = @(r,theta,delta) R0+delta+r.*cos(theta);
    D = @(r,theta,k,kprime,deltaprime) kprime.*r.*sin(theta).^2+k.*deltaprime.*cos(theta)+k;
    BdotGradPhi =  @(r,theta,delta) B0*R0./(R(r,theta,delta).^2);
    BdotGradTheta = @(r,theta,k,kprime,deltaprime,delta,qbar) B0.*k./(qbar.*R(r,theta,delta).*D(r,theta,k,kprime,deltaprime));
    
    integrand = @(theta) BdotGradPhi(r,theta,delta)./BdotGradTheta(r,theta,k,kprime,deltaprime,delta,qbar);
    t=zeros(1,npoints);
    
    for i=1:npoints
        t(i) = (1/q)*midpoint_composite_quadrature(integrand,0,th(i),i-1) ; 
    end
 
return