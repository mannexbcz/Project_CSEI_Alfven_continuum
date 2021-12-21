function qbar = qbar_ordre1(q,r,k,kprime)
% This function implements the numerical computation of qbar in the first
% order approximation, using equation (58)
% Inputs:
%       q : safety factor, handle function
%       r : radius, double
%       k, kprime : coefficients for the elongation, double
% Outputs:
%       qbar, double

qbar = (2.*k.*q)./(kprime.*r+2.*k);
return