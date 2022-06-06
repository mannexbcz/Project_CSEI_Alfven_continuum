function [M,N]=build_MN_shifted(epsilon,delta_p,n,m_min, m_max, q,F,size)
% This function build the matrices M and N in the case of a circular, 
% and shifted equilibriuum. The form of the matrices are given by equations
% (39) ans (40) of the report.
% Inputs:
%   epsilon : large aspect ratio, double
%   delta_p : derivative of the shift, double
%   n : toroidal mode number, int
%   m_min, m_max : range of values for the poloidal mode number m
%   q : safety factor, handle function
%   F : toroidal flux current, double
%   size : required rank of the matrix M, int
% Output: 
%   M : matrix M, of size (size,size)
%   N : matrix N, of size (size,size)

    m=[m_min:1:m_max];

    M0=ones(1,size);
    M1=(2*epsilon-delta_p)*ones(1,size-1);
    M=diag(M0)+diag(M1,-1)+diag(M1,1);
    
    N0=ones(1,size);
    N1=-delta_p*ones(1,size-1);
    N=diag(N0)+diag(N1,-1)+diag(N1,1);
    
    for i=1:size
        for j=1:size
            N(i,j)=N(i,j)*F^2*(n+m(i)./q).*(n+m(j)./q);
        end
    end

end