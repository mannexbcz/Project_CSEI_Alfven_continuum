% -------------------------------------------------------------------------
% This script compares the Fourier coefficients obtained from the
% equilibrium coefficients directly given by CHEASE and the ones
% parametrized (using Delta,Kappa,delta etc provided by CHEASE)
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;
%%

load('coeffgeom.mat')
load('coeffchease.mat')
load('coeffchease_redone.mat')

coeffsMchease = real(coeffsMchease);
coeffsNchease = real(coeffsNchease);
coeffsMgeom = real(coeffsMgeom);
coeffsNgeom = real(coeffsNgeom);
coeffsMchease_redone = real(coeffsMchease_redone);
coeffsNchease_redone = real(coeffsNchease_redone);

filename = 'ogyropsi_PT_redone_nmeshpol_0.h5';
params = read_chease(filename);
rchease_redone = params.r;

filename = 'ogyropsi_positive_triangularity.h5';
params = read_chease(filename);
rchease = params.r;
%% comparaison chease/geom

save=1; %Set 1 to save the figure, 0 otherwise
namefig = 'Comp_Fourier_M_CHEASE_GEOM';
path='C:\Users\manon\Desktop\projet CSE II\figures\';

figure
plot(rchease, coeffsMchease(:,2)./coeffsMchease(:,1), 'k-')
hold on 
plot(rchease, coeffsMchease(:,3:end)./coeffsMchease(:,1), 'k-', 'HandleVisibility','off')
hold on
plot(rchease,coeffsMgeom(:,2)./coeffsMgeom(:,1), 'k--')
hold on
plot(rchease,coeffsMgeom(:,3:end)./coeffsMgeom(:,1), 'k--','HandleVisibility','off')
xlabel('$r$')
ylabel('Fourier coefficients of M')
legend('CHEASE','Geometric')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

save=1; %Set 1 to save the figure, 0 otherwise
namefig = 'Comp_Fourier_N_CHEASE_GEOM';
path='C:\Users\manon\Desktop\projet CSE II\figures\';

figure
plot(rchease, coeffsNchease(:,2)./coeffsNchease(:,1), 'k-')
hold on 
plot(rchease, coeffsNchease(:,3:end)./coeffsNchease(:,1), 'k-', 'HandleVisibility','off')
hold on
plot(rchease,coeffsNgeom(:,2)./coeffsNgeom(:,1), 'k--')
hold on
plot(rchease,coeffsNgeom(:,3:end)./coeffsNgeom(:,1), 'k--','HandleVisibility','off')
xlabel('$r$')
ylabel('Fourier coefficients of N')
legend('CHEASE','Geometric')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end


%% Mean value of coefficients

meanMchease = mean(coeffsMchease);
meanNchease = mean(coeffsNchease);
meanMgeom = mean(coeffsMgeom);
meanNgeom = mean(coeffsNgeom);

x = [1:7];
figure
bar(x,[meanMchease;meanNchease])
xlabel('Fourier coefficient')
ylabel('Mean amplitude')
legend('M','N')

%% comparaison chease/redone

save=1; %Set 1 to save the figure, 0 otherwise
namefig = 'Comp_Fourier_M_CHEASE_redone';
path='C:\Users\manon\Desktop\projet CSE II\figures\';

figure
plot(rchease, coeffsMchease(:,2)./coeffsMchease(:,1), 'k-')
hold on 
plot(rchease, coeffsMchease(:,3:end)./coeffsMchease(:,1), 'k-', 'HandleVisibility','off')
hold on
plot(rchease_redone,coeffsMchease_redone(:,2)./coeffsMchease_redone(:,1), 'k--')
hold on
plot(rchease_redone,coeffsMchease_redone(:,3:end)./coeffsMchease_redone(:,1), 'k--','HandleVisibility','off')
xlabel('$r$')
ylabel('Fourier coefficients of M')
legend('non-equidistant $\theta^*$','equidistant $\theta^*$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

save=1; %Set 1 to save the figure, 0 otherwise
namefig = 'Comp_Fourier_N_CHEASE_redone';
path='C:\Users\manon\Desktop\projet CSE II\figures\';

figure
plot(rchease, coeffsNchease(:,2)./coeffsNchease(:,1), 'k-')
hold on
plot(rchease, coeffsNchease(:,3:end)./coeffsNchease(:,1), 'k-', 'HandleVisibility','off')
hold on
plot(rchease_redone,coeffsNchease_redone(:,2)./coeffsNchease_redone(:,1), 'k--')
hold on
plot(rchease_redone,coeffsNchease_redone(:,3:end)./coeffsNchease_redone(:,1), 'k--','HandleVisibility','off')
xlabel('$r$')
ylabel('Fourier coefficients of N')
legend('non-equidistant $\theta^*$','equidistant $\theta^*$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end