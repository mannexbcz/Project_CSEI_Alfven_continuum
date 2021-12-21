% -------------------------------------------------------------------------
% This script verifies the convergence of the computed Fourier coefficients
% with respect to the number of mesh points used for their integration.
% The case of an elongated and shifted equilibriuum, in the first order 
% approximation is used.
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
nr = 1000;  % number of magnetic surfaces considered
npoints=100;% number of mesh points used for the integration
q = @(r) 1+2*(r/a).^2; % safety factor

% Elongation
k_a=1.2;
kprime=0.1;
kfct = @(r) kprime.*(r-a) + k_a;

%Shift
deltaprime= -0.1;
deltafct = @(r)  -deltaprime*a + deltaprime *r;

qbarr = @(r,k,kprime) (2.*k.*q(r))./(kprime.*r+2.*k);

%%
% Useful functions
thetastar = @(r,theta,k,kprime,deltaprime,epsilon) theta+(deltaprime-epsilon).*sin(theta)+(kprime.*r./(4.*k))*sin(2.*theta);
dThetastardTheta = @(r,theta,k,kprime,deltaprime,delta,epsilon,qbar) (kprime.*r./k).*sin(theta).^2+(deltaprime-epsilon).*cos(theta)+1;

% Equilibriuum coefficients
Mtheta = @(r,theta,k,kprime,delta,deltaprime,epsilon,qbar) (B0^4*R0^4)./(qbar.^2).*(k.^2.*cos(theta).^2+sin(theta).^2).*(1-(2*kprime*R0/k).*epsilon.*sin(theta).^2-2.*deltaprime.*cos(theta)).*(1+4.*epsilon.*cos(theta));
Ntheta = @(r,theta,k,kprime,delta,deltaprime,epsilon,qbar) (B0^2*R0^2)./(qbar.^2).*(k.^2.*cos(theta).^2+sin(theta).^2).*(1-(2*kprime*R0/k).*epsilon.*sin(theta).^2-2.*deltaprime.*cos(theta));

%% Reference

 coeffs_M_REF=[];
 coeffs_N_REF=[];
 
r=0.05;
epsilon=r/R0;
delta = deltafct(r);
k=kfct(r);
qbar=qbarr(r,k,kprime);
nref=10^6;

    for m = -10:10
        fun_theta1 = @(theta) exp(-1i.*m.*thetastar(r,theta,k,kprime,deltaprime,epsilon)).*Mtheta(r,thetastar(r,theta,k,kprime,deltaprime,epsilon),k,kprime,delta,deltaprime,epsilon,qbar).*dThetastardTheta(r,theta,k,kprime,deltaprime,delta,epsilon,qbar);
        c1=(1/2*pi)*midpoint_composite_quadrature(fun_theta1, 0, 2*pi, nref);
        coeffs_M_REF=[coeffs_M_REF,c1];
        
        fun_theta2 = @(theta) exp(-1i.*m.*thetastar(r,theta,k,kprime,deltaprime,epsilon)).*Ntheta(r,thetastar(r,theta,k,kprime,deltaprime,epsilon),k,kprime,delta,deltaprime,epsilon,qbar).*dThetastardTheta(r,theta,k,kprime,deltaprime,delta,epsilon,qbar);
        c2=(1/2*pi)*midpoint_composite_quadrature(fun_theta2, 0, 2*pi, nref);
        coeffs_N_REF=[coeffs_N_REF,c2];
        
    end

%% 

ns=round(logspace(0,5,20));
errorM=zeros(1,length(ns));
errorN=zeros(1,length(ns));

for s=1:length(ns)
    
npoints=ns(s);

r=0.05;

    [coeffs_M] = get_fourier_coeff_theta_prime(Mtheta,r,k,kprime,delta,deltaprime,epsilon,qbar,thetastar,dThetastardTheta,-10,10,npoints);
    [coeffs_N] = get_fourier_coeff_theta_prime(Ntheta,r,k,kprime,delta,deltaprime,epsilon,qbar,thetastar,dThetastardTheta,-10,10,npoints);
    
    errorM(s)=max(abs((coeffs_M-coeffs_M_REF)/coeffs_M_REF));
    errorN(s)=max(abs((coeffs_N-coeffs_N_REF)/coeffs_N_REF));

end

%% FIGURE

save=1;%Set 1 to save the figure, 0 otherwise
namefig=['Convergence_coeff_fourier_first_order'];
path='C:\Users\manon\Desktop\projet CSE I\figures\elongation\';

figure
loglog(ns,errorM,'b+')
hold on
loglog(ns,errorN,'r+')
grid on
xticks([10^0 10^1 10^2 10^3 10^4 10^5])
xticklabels({'$10^0$','$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'})
xlabel('Number of mesh points')
ylabel('Relative Error')
legend('$a_\mathrm{eq}$','$b_\mathrm{eq}$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end