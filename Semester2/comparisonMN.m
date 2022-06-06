% -------------------------------------------------------------------------
% This script compares the coefficients M and N obtained with and without
% an equidistant mesh in Thetastar
% -------------------------------------------------------------------------

%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;
%%
load('geom.mat')
load('chease.mat')
load('chease_redone.mat')

filename = 'ogyropsi_PT_redone_nmeshpol_0.h5';
params = read_chease(filename);
r_redone = params.r;
thetastar_chease_redone = params.thetastar;

filename = 'ogyropsi_positive_triangularity.h5';
params = read_chease(filename);
r = params.r;
thetastar_chease = params.thetastar;

thetastar_geom = zeros(3,length(thetastar_chease));
npoints = length(thetastar_chease);
theta = linspace(0,2*pi,npoints);
j=0;
for i=[175,206,223]
    j=j+1;
    epsilon=params.r(i)/params.R0;
    delta = params.Delta(i);
    k = params.kappa(i);
    d = params.delta(i);
    kprim = params.kprime(i);
    deltaprim = params.deltaprime(i);
    dprim = params.dprime(i);
    qbar = qbar_triang(params.q(i),params.r(i),k,kprim,delta,deltaprim,d,dprim,params.R0,npoints);
    
    thetastar_geom(j,:) = theta_star_triang(r(i),theta,kprim,k,delta,deltaprim,d,dprim,params.q(i),qbar,params.B0,params.R0,npoints); 
end

%% M
save=0; %Set 1 to save the figure, 0 otherwise
namefig = 'Mq15';
path='C:\Users\manon\Desktop\projet CSE II\figures\';

figure
plot(thetastar_chease, Mchease(175,:))
hold on
plot(thetastar_geom(1,:), Mgeom(175,:))
xlabel('$\theta^*$')
ylabel('$M$')
title('$q(r)=1.5$')
legend('CHEASE','Geometric')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

%% N
save=0; %Set 1 to save the figure, 0 otherwise
namefig = 'Nq15';
path='C:\Users\manon\Desktop\projet CSE II\figures\';

figure
plot(thetastar_chease, Nchease(175,:))
hold on
plot(thetastar_geom(1,:), Ngeom(175,:))
xlabel('$\theta^*$')
ylabel('$N$')
title('$q(r)=1.5$')
legend('CHEASE','Geometric')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

%% M
save=1; %Set 1 to save the figure, 0 otherwise
namefig = 'Mq15_equid_non_equid';
path='C:\Users\manon\Desktop\projet CSE II\figures\';

figure
plot(thetastar_chease, Mchease(175,:))
hold on
plot(thetastar_chease_redone, Mchease_redone(175,:))
xlabel('$\theta^*$')
ylabel('$M$')
title('$q(r)=1.5$')
legend('non-equidistant $\theta^*$','equidistant $\theta^*$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

%% N
save=1; %Set 1 to save the figure, 0 otherwise
namefig = 'Nq15_equid_non_equid';
path='C:\Users\manon\Desktop\projet CSE II\figures\';

figure
plot(thetastar_chease, Nchease(175,:))
hold on
plot(thetastar_chease_redone, Nchease_redone(175,:))
xlabel('$\theta^*$')
ylabel('$N$')
title('$q(r)=1.5$')
legend('non-equidistant $\theta^*$','equidistant $\theta^*$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end