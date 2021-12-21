function w = eigenmodes_triang(size,band,kfct,k_a,kprime,deltafct,deltaprime,dfct,dprime,q,R0,B0,a,n,nr,npoints)
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

r=linspace(0,a,nr);
w=zeros(size,nr);
w2=zeros(size,nr);

for i=1:nr
    epsilon=r(i)/R0;
    delta = deltafct(r(i));
    k = kfct(r(i));
    d = dfct(r(i));
    
    qbar = qbar_triang(q(r(i)),r(i),k,kprime,delta,deltaprime,d,dprime,R0,npoints);

    [M,N] = matrices_triang(r(i),epsilon,k,kprime,delta,deltaprime,d,dprime,q(r(i)),qbar,R0,B0,n,size,band,npoints);

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