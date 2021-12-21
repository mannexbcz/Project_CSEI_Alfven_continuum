function B=build_N(n,m_min, m_max, q,F,size)
% This function build the matrix N in the case of a circular, concentric
% and large aspect ratio equilibriuum. Its form is given by equation (29)
% of the report.
% Inputs:
%   n : toroidal mode number, int
%   m_min, m_max : range of values for the poloidal mode number m
%   size : required rank of the matrix M, int 
%   q : safety factor, handle function
%   F : toroidal flux current, double
% Output: 
%   N : matrix N, of size (size,size)

    N0 = zeros(1,size);
    m = [m_min:1:m_max];
    for i=1:size
        N0(i)=F^2.*(n+m(i)/q)^2;
    end
    
    B=diag(N0);   
end