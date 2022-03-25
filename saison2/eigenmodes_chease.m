function w = eigenmodes_chease(size,band,gradPsi2,B,B0,R,R0,F,q,thetastar,a,n,nr)
% Computes the Alfvén eigenmodes by solving the matrix form of the Alfvén
% continuum equation, for the elongated and shifted case. 
% Inputs:
%   size : required rank of the matrices, int
%   band : required number of bands of the matrices, int
%   r : radius, double
%   kfct,k_a, kprime : parametrisation of the elongation, handle
%                       function,double,double
%   deltafct, deltaprime :parametrisation of the shift, handle function, double
%   q : safety factor, handle function
%   n : toroidal mode number, int
%   R0 : major radius, double
%   B0 : magnetic fiels, double
%   nr : number of magnetic surfaces
%   npoints : number of mesh points for the numerical integration, int
% Output: 
%   w : Alfvén eigenmodes, matrix of size (size,nr)

% r=transpose(linspace(0,a,nr));
w=zeros(size,nr);
w2=zeros(size,nr);

M_all = gradPsi2.*(B0./B).^2.*(R./R0).^2;
N_all = gradPsi2./(B.^2.*R.^2);

% Mchease = M_all(1,:);
% Nchease = N_all(1,:);
% save('chease.mat','Mchease','Nchease')

for i=1:nr
    
    [M,N] = matrices_chease(M_all(i,:),N_all(i,:),thetastar,q(i),F(i),n,size,band);

    w2(:,i)=eig(N,M); 

    for j=1:size 
        if w2(j,i)<0
            w2(j,i)=nan;
        end
    end
    w2(:,i)=sort(w2(:,i));
    w(:,i)=sqrt(w2(:,i));
end

return