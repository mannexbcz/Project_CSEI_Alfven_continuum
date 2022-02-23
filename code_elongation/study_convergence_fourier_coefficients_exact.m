% -------------------------------------------------------------------------
% This script verifies the convergence of the computed Fourier coefficients
% with respect to the number of mesh points used for their integration.
% The case of an elongated and shifted equilibriuum is used.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;
%% Parameters

size=5;     % rank of the matrices M and N
band=5;     % bands of the matrices M and N
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
B0=1;       % magnetic field
F = R0*B0;  % toroidal flux current
nr = 100;  % number of magnetic surfaces considered
npoints=100;% number of mesh points used for the integration
qfct = @(r) 1+2*(r/a).^2; % safety factor

% Elongation
k_a=1.2;
kprime=0.1;
kfct = @(r) kprime.*(r-a) + k_a;

%Shift
deltaprime= -0.1;
deltafct = @(r)  -deltaprime*a + deltaprime *r;

% Useful fonctions
GradPsi2 = @(r,theta,k,kprime,deltaprime,qbar) B0^2.*r.^2.*k.^2.*((k.^2.*cos(theta).^2+sin(theta).^2))./(qbar.^2.*(k.*deltaprime.*cos(theta)+kprime.*r.*sin(theta).^2+k).^2); %B0^2.*r.^2.*k.^2.*
R = @(r,theta,delta) R0+delta+r.*cos(theta);
B2 = @(r,theta,k,kprime,delta,deltaprime,qbar) (1./R(r,theta,delta).^2).*(R0^2*B0^2+GradPsi2(r,theta,k,kprime,deltaprime,qbar));
thetastar = @(r,theta,k,kprime,deltaprime,epsilon) (theta+(deltaprime-epsilon).*sin(theta)+(kprime.*r./(4.*k))*sin(2.*theta));
D = @(r,theta,k,kprime,deltaprime) kprime.*r.*sin(theta).^2+k.*deltaprime.*cos(theta)+k;
BdotGradPhi =  @(r,theta,delta) B0*R0./(R(r,theta,delta).^2);
BdotGradTheta = @(r,theta,k,kprime,deltaprime,delta,qbar) B0.*k./(qbar.*R(r,theta,delta).*D(r,theta,k,kprime,deltaprime));
dThetastardTheta = @(r,theta,k,kprime,deltaprime,delta,epsilon,qbar) BdotGradPhi(r,theta,delta)./BdotGradTheta(r,theta,k,kprime,deltaprime,delta,qbar);

% Equilibriuum coefficients
Mtheta = @(r,theta,k,kprime,delta,deltaprime,epsilon,qbar) (B0^2/R0^2).*((R(r,theta,delta).^2)./B2(r,theta,k,kprime,delta,deltaprime,qbar)).*GradPsi2(r,theta,k,kprime,deltaprime,qbar) ;
Ntheta = @(r,theta,k,kprime,delta,deltaprime,epsilon,qbar) GradPsi2(r,theta,k,kprime,deltaprime,qbar)./(B2(r,theta,k,kprime,delta,deltaprime,qbar).*(R(r,theta,delta).^2)) ;


%% Reference

r=0.05;
delta=deltafct(r);
k=kfct(r);
epsilon=r/R0;
q = qfct(r);

nref=10^4;
mmin=0;
mmax=10;

qbar = qbar_exact(q,r,epsilon,k,kprime,delta,deltaprime,R0,nref);

th=linspace(0,2*pi,nref);
thetastar=theta_star(r,th,kprime,k,delta,deltaprime,q,qbar,B0,R0,nref);
thstarmid=0.5*(thetastar(1:nref-1)+thetastar(2:nref));

[coeffs_M_REF] = get_fourier_coeff_theta_prime_exact(Mtheta,r,k,kprime,delta,deltaprime,epsilon,q,qbar,thstarmid,dThetastardTheta,mmin,mmax,nref,R0,B0);
[coeffs_N_REF] = get_fourier_coeff_theta_prime_exact(Ntheta,r,k,kprime,delta,deltaprime,epsilon,q,qbar,thstarmid,dThetastardTheta,mmin,mmax,nref,R0,B0);
    

%%

ns=round(logspace(0,4,20));
errorM=zeros(1,length(ns));
errorN=zeros(1,length(ns));

for s=1:length(ns)
    
npoints=ns(s);

r=0.05;
qbar = qbar_exact(q,r,epsilon,k,kprime,delta,deltaprime,R0,npoints);

    th=linspace(0,2*pi,npoints);
    thetastar=theta_star(r,th,kprime,k,delta,deltaprime,q,qbar,B0,R0,npoints);
    thstarmid=0.5*(thetastar(1:npoints-1)+thetastar(2:npoints));

    [coeffs_M] = get_fourier_coeff_theta_prime_exact(Mtheta,r,k,kprime,delta,deltaprime,epsilon,q,qbar,thstarmid,dThetastardTheta,mmin,mmax,npoints,R0,B0);
    [coeffs_N] = get_fourier_coeff_theta_prime_exact(Ntheta,r,k,kprime,delta,deltaprime,epsilon,q,qbar,thstarmid,dThetastardTheta,mmin,mmax,npoints,R0,B0);
     
    errorM(s)=max(abs((coeffs_M-coeffs_M_REF)/coeffs_M_REF));
    errorN(s)=max(abs((coeffs_N-coeffs_N_REF)/coeffs_N_REF));

end

%% FIGURE

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['Convergence_coeff_fourier_exact'];
path='C:\Users\manon\Desktop\projet CSE I\figures\elongation\';

Const1=polyfit(log(ns(5:end-1)),log(errorM(5:end-1)),1);
m = Const1(1);
k= Const1(2);
JBL = ns(1:end-1).^m.*exp(k);

Const2=polyfit(log(ns(5:end-1)),log(errorN(5:end-1)),1);
m2 = Const2(1);
k2= Const2(2);
JBL2 = ns(1:end-1).^m2.*exp(k2);

figure
loglog(ns,errorM,'b+')
hold on
loglog(ns(1:end-1),JBL,'b--')
hold on
loglog(ns,errorN,'r+')
hold on
loglog(ns(1:end-1),JBL2,'r--')
grid on
xticks([10^0 10^1 10^2 10^3 10^4 10^5])
xticklabels({'$10^0$','$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'})
xlabel('Number of mesh points')
ylabel('Relative Error')
strlegend1=['$y=$',num2str(m,3),'$x+b$'];
strlegend2=['$y=$',num2str(m2,3),'$x+b$'];
legend('$a_\mathrm{eq}$',strlegend1,'$b_\mathrm{eq}$',strlegend2)

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end