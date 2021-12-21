function qbar = qbar_exact(q,r,epsilon,k,kprime,delta,deltaprime,R0,npoints)
% This function implements the numerical computation of qbar according to
% equation (21)
% Inputs:
%       q : safety factor, handle function
%       r : radius, double
%       epsilon : inverse aspect ratio, double,
%       k, kprime : coefficients for the elongation, double
%       delta, deltaprime : coefficients for the shift, double
%       R0 : major radius, double
%       npoints : number of mesh points for the numerical integration, int
% Outputs:
%       qbar, double

    integrand_q = @(theta)  ((kprime.*r./k).*sin(theta).^2+deltaprime.*cos(theta)+1)./(1+delta./R0+epsilon.*cos(theta));

    qbar = q*2*pi/midpoint_composite_quadrature(integrand_q, 0, 2*pi, npoints);
end