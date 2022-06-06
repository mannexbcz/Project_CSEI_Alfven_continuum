function [M,N] = matrices_chease(Mi,Ni,thetastar,q,F,n,size,band,i)
% This function build the matrices M and N for the solution obtained from
% the equilibrium coefficients provided by CHEASE
% Inputs:
%   Mi,Ni: vectors, equilibrium coefficients at a given magnetic surface
%   q : safety factor, doubles
%   thetastar : vector of doubles
%   n : toroidal mode number, int
%   F : toroidal current flux function, double
%   size : required rank of the matrices, int
%   band : required number of bands of the matrices, int
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
    
        savefourier = 0;
if savefourier
    load('coeffchease.mat')
    coeffsMchease(i,:) = coeffs_M;
    coeffsNchease(i,:) = coeffs_N;
    save('coeffchease.mat','coeffsMchease','coeffsNchease')
end
    
   %Construction of the matrices
    [M,N]=build_MN_bands_Fourier(coeffs_M, coeffs_N, m_min,m_max, n, q, F, size,band);
return


