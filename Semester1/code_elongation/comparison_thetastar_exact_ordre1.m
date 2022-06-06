% -------------------------------------------------------------------------
% This scripts compares the values of theta star computed using the exact
% integration and the first order approximation.
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
r=0.05;
q = 1+2*(r/a).^2; % safety factor

% Elongation
k_a=1.1;
kprime=0;
kfct = kprime.*(r-a) + k_a;

%Shift
deltaprime= -0.1;
deltafct =  -deltaprime*a + deltaprime *r;

epsilon=r/R0;

%%

integrand_q = @(theta)  ((kprime.*r./k).*sin(theta).^2+deltaprime.*cos(theta)+1)./(1+delta./R0+epsilon.*cos(theta));
qbar = q*2*pi/midpoint_composite_quadrature(integrand_q, 0, 2*pi, npoints);

th = linspace(0,2*pi,npoints);
thetastar = @(r,theta,k,kprime,deltaprime,epsilon) (theta+(deltaprime-epsilon).*sin(theta)+(kprime.*r./(4.*k))*sin(2.*theta));

texact = theta_star(r,th,kprime,k,delta,deltaprime,q,qbar,B0,R0,npoints);
tordre1 = thetastar(r,th,k,kprime,deltaprime,epsilon);

%% Figure : theta star as a function of theta, for the two implementations

figure
plot(th,texact)
hold on
plot(th,tordre1)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'})
yticks([0 pi/2 pi 3*pi/2 2*pi])
yticklabels({'0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'})
xlabel('$\theta$')
ylabel('$\theta^*$')
legend('exact','$1^\mathrm{st}$ order')
grid on

%% Figure : error between the two implementations of theta star

err=abs((texact-tordre1)./texact);

figure
semilogy(th,err)
grid on
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'})
xlabel('$\theta$')
ylabel('Relative error')

%% Error as a function of npoints
number = 10;
ns = round(logspace(1,5,number));
errors = zeros(1,number);

for i = 1: number
    npoints=ns(i);
    qbar = q*2*pi/midpoint_composite_quadrature(integrand_q, 0, 2*pi, npoints);
    th = linspace(0,2*pi,npoints);
    thetastar = @(r,theta,k,kprime,deltaprime,epsilon) (theta+(deltaprime-epsilon).*sin(theta)+(kprime.*r./(4.*k))*sin(2.*theta));

    texact = theta_star(r,th,kprime,k,delta,deltaprime,q,qbar,B0,R0,npoints);
    tordre1 = thetastar(r,th,k,kprime,deltaprime,epsilon);
    errors(i) = max(abs((texact-tordre1)./tordre1));
end

figure
loglog(ns,errors,'+ ')
grid on
xlabel('$n$')
ylabel('Relative error')
