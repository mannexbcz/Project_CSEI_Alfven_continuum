function [M,N] = matrices_chease(Mi,Ni,thetastar,q,F,n,size,band)
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

    m_min=round((-n*q-size/2+1/2));
    m_max=m_min+size+1;

    m_min_coeff=0;
    m_max_coeff=(band+1)/2;
    
   % Fourier coefficients(fun,thetastar,mmin,mmax)
    [coeffs_M] = get_fourier_coeff_chease(Mi,thetastar,m_min_coeff,m_max_coeff);
    [coeffs_N] = get_fourier_coeff_chease(Ni,thetastar,m_min_coeff,m_max_coeff);
    
   %Construction of the matrices
    [M,N]=build_MN_bands_Fourier(coeffs_M, coeffs_N, m_min,m_max, n, q, F, size,band);
return


