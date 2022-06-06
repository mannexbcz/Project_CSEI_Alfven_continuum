function w = eigenmodes_chease_geom(size,band,params,n,nr,npoints)
% Computes the Alfvén eigenmodes by solving the matrix form of the Alfvén
% continuum equation, for geometrically parametrized equilibria, with
% coefficients provided by CHEASE
% Inputs:
%   size : required rank of the matrices, int
%   band : required number of bands of the matrices, int
%   params : parameters read from CHEASE (.h5 files)
%   n : toroidal mode number, int
%   nr : number of magnetic surfaces
%   npoints : number of mesh points for the numerical integration, int
% Output: 
%   w : Alfvén eigenmodes, matrix of size (size,nr)

% r=transpose(linspace(0,a,nr));
w=zeros(size,nr);
w2=zeros(size,nr);


savecoeff=0;
if savecoeff
Mgeom = zeros(length(params.r),length(params.thetastar));
Ngeom = zeros(length(params.r),length(params.thetastar));
save('geom.mat','Mgeom','Ngeom')
end

savefourier=0;
if savefourier
coeffsMgeom = [];
coeffsNgeom = [];
save('coeffgeom.mat','coeffsMgeom','coeffsNgeom')
end

for i=1:nr
    epsilon=params.r(i)/params.R0;
    delta = params.Delta(i);
        k = params.kappa(i);
        d = params.delta(i);
        kprim = params.kprime(i);
        deltaprim = params.deltaprime(i); 
        dprim = params.dprime(i);
        qbar = qbar_triang(params.q(i),params.r(i),k,kprim,delta,deltaprim,d,dprim,params.R0,npoints);
    
    [M,N] = matrices_triang(params.r(i),epsilon,k,kprim,delta,deltaprim,d,dprim,params.q(i),qbar,params.R0,params.B0,n,size,band,npoints,i);

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