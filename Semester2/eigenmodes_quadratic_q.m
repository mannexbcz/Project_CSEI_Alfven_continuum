function w = eigenmodes_quadratic_q(size,band,params,n,nr,npoints)
% Computes the Alfvén eigenmodes by solving the matrix form of the Alfvén
% continuum equation, for the quadratic approximation of $q$ 
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

% Computation of deltaq, r1 
Deltaq = 1-params.q(1) ;
disp(Deltaq)
r1 = params.a *(1-(1-params.q(end)) /Deltaq)^(-1/2);
% Computation of s(r), delta(r), Kappa(r)

s = (2*Deltaq./params.q).*(params.r./r1).^2;
Delta0 = params.Delta(1); 
ka = params.kappa(end); sa = s(end) ; da = params.delta(end); qa = params.q(end);
kappafct = ka - ((ka^2-1).*(qa*sa-params.q.*s))./(12+(1+ka)*qa*sa-(ka-1).*params.q.*s);
deltafct = da.*(params.r./params.a).*(4+params.q.*s)./(4+qa*sa);
kappaprimefct = derive(kappafct,params.r);
deltaprimefct = derive(deltafct,params.r);
Deltaprimefct = Deltaprime_quadratic(params.r,params.R0,params.F,params.q,params.pprime);
Deltafct = Delta_quadratic(Deltaprimefct,params.r,Delta0);
q_quad =  1-Deltaq.*(1-(params.r./r1).^2);
figure
plot(params.r,Deltaprimefct)
xlabel('$r$')
ylabel('$\Delta^\prime$')

figure
plot(params.r,params.Delta)
hold on
plot(params.r,Deltafct)
legend('CHEASE','Quadratic approximation of $q$')
xlabel('$r$')
ylabel('$\Delta$')

figure
plot(params.r,params.kappa)
hold on
plot(params.r,kappafct)
legend('CHEASE','Quadratic approximation of $q$')
xlabel('$r$')
ylabel('$\kappa$')

figure
plot(params.r,params.delta)
hold on
plot(params.r, deltafct)
legend('CHEASE','Quadratic approximation of $q$')
xlabel('$r$')
ylabel('$\delta$')

figure
plot(params.r./params.a,params.q)
hold on
plot(params.r./params.a, 1-Deltaq.*(1-(params.r./r1).^2))
legend('CHEASE','Quadratic approximation')
xlabel('$r$')
ylabel('$q$')

for i=1:nr
    epsilon=params.r(i)/params.R0;
    delta = Deltafct(i);
        k = kappafct(i);
        d = deltafct(i);
        kprim = kappaprimefct(i);
        deltaprim = Deltaprimefct(i); 
        dprim =deltaprimefct(i);
        qbar = qbar_triang(q_quad(i),params.r(i),k,kprim,delta,deltaprim,d,dprim,params.R0,npoints);
        
        %disp(((ka^2-2).*(qa*sa-params.q.*s)))
%         disp('test2')
%         disp(delta)
%         disp(deltaprim)
    [M,N] = matrices_triang(params.r(i),epsilon,k,kprim,delta,deltaprim,d,dprim,q_quad(i),qbar,params.R0,params.B0,n,size,band,npoints,i);

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