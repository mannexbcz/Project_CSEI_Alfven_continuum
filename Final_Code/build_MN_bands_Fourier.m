function [Mmat,Nmat]=build_MN_bands_Fourier(coeffsM, coeffsN, m_min, m_max, n, q, F, size, band)
% This function build the matrices M and N from the Fourier coefficients of 
% coefficients M and N, applying equations (17-19)
% Inputs:
%   coeffsM : Fourier coefficients of equilibrium coefficient M
%   coeffsN : Fourier coefficients of equilibrium coefficient N
%   n : toroidal mode number, int
%   m_min, m_max : range of values for the poloidal mode number m
%   q : safety factor
%   F : toroidal flux current, double
%   size : required rank of the matrices, int
%   band : required number of bands of the matrices, int
% Output: 
%   Mmat : matrix M, of size (size,size)
%   Nmat : matrix N, of size (size,size)
    
    Mmat=zeros(size,size);
    Nmat=zeros(size,size);
    m=[m_min:1:m_max];

    idx_max=(band-1)/2;
    for i=1:size
        for jeq=0:idx_max
            j = i+jeq;
            if j<=size
                Mmat(i,j)=coeffsM(jeq+1);
                Nmat(i,j)=coeffsN(jeq+1)*F^2*(n+m(i)./q).*(n+m(j)./q);
                Mmat(j,i)=Mmat(i,j);
                Nmat(j,i)=Nmat(i,j);
            end
        end
    end

end