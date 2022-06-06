function w = eigenmodes_chease(size,band,gradPsi2,B,B0,R,R0,F,q,thetastar,a,n,nr)
% Computes the Alfvén eigenmodes by solving the matrix form of the Alfvén
% continuum equation, for the coefficients provided by CHEASE. 
% Inputs:
%   size : required rank of the matrices, int
%   band : required number of bands of the matrices, int
%   gradPsi2 : provided by CHEASE, vector of doubles
%   q : safety factor, vector of doubles
%   n : toroidal mode number, int
%   R0 : major radius, double
%   B0 : magnetic fields, double
%   R : vector of doubles
%   B : vector of doubles
%   nr : number of magnetic surfaces
%   F : toroidal current flux function, vector of doubles
%   a : double
% Output: 
%   w : Alfvén eigenmodes, matrix of size (size,nr)

w=zeros(size,nr);
w2=zeros(size,nr);

M_all = gradPsi2.*(B0./B).^2.*(R./R0).^2;
N_all = gradPsi2./(B.^2.*R.^2);

Mchease = M_all(:,:);
Nchease = N_all(:,:);
save('chease.mat','Mchease','Nchease')

savefourier=0;
if savefourier
coeffsMchease = [];
coeffsNchease = [];
save('coeffchease.mat','coeffsMchease','coeffsNchease')
end

for i=1:nr
    
    [M,N] = matrices_chease(M_all(i,:),N_all(i,:),thetastar,q(i),F(i),n,size,band,i);

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