function M=build_M(epsilon, size)
% This function build the matrix M in the case of a circular, concentric
% and large aspect ratio equilibriuum. Its form is given by equation (28)
% of the report.
% Inputs:
%   epsilon : large aspect ratio, double
%   size : required rank of the matrix M, int
% Output: 
%   M : matrix M, of size (size,size)
    
    M0 = ones(1,size);
    M1 = 2*epsilon*ones(1,size-1);
    M = diag(M0)+diag(M1,-1)+diag(M1,1);
    
end