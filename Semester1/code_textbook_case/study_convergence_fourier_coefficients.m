% -------------------------------------------------------------------------
% This script verifies the convergence of the computed Fourier coefficients
% with respect to the number of mesh points used for their integration.
% The textbook case (eq.(78)) is used.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;
%% Parameters

size=5;             % rank of the matrices M and N
N_modes = 3;        % Number of modes 
band=2*N_modes-1;   % number of bands of the matrices M and N
a=0.1;              % minor radius
n=1;                % toroidal mode number
R0=1;               % major radius
B0=1;               % magnetic field
F = R0*B0;          % toroidal flux current
nr = 100;          % number of magnetic surfaces considered
q = @(r) 1+2*(r/a).^2; % safety factor

%% Equilibrium coefficients

M = @(epsilon,epsilon2,epsilon3,theta) 1+4.*epsilon.*cos(theta+epsilon2.*sin(theta));
N = @(epsilon,epsilon2,epsilon3,theta) 1+epsilon3.*cos(theta+epsilon2.*sin(theta)).^2;

epsilon = 0.1;
epsilon2 = 0.2;
epsilon3 = 0.3;

%% Reference

 nref=10^6;     % number of mesh points for the numerical integration
 
[coeffs_M_REF] = get_fourier_coeff_pt_milieu(M,epsilon,epsilon2,epsilon3,-10,10,nref);
[coeffs_N_REF] = get_fourier_coeff_pt_milieu(N,epsilon,epsilon2,epsilon3,-10,10,nref);

%% Loop over the number of mesh points

ns=round(logspace(0,5,20));  % number of mesh points for the numerical integration
errorM=zeros(1,length(ns));
errorN=zeros(1,length(ns));

for s=1:length(ns)
    
npoints=ns(s);

[coeffs_M] = get_fourier_coeff_pt_milieu(M,epsilon,epsilon2,epsilon3,-10,10,npoints);
[coeffs_N] = get_fourier_coeff_pt_milieu(N,epsilon,epsilon2,epsilon3,-10,10,npoints);

errorM(s)=max(abs((coeffs_M-coeffs_M_REF)/coeffs_M_REF));
errorN(s)=max(abs((coeffs_N-coeffs_N_REF)/coeffs_N_REF));
end

%%
save=0; %Set 1 to save the figure, 0 otherwise
namefig=['cas_ecole_convergence_coeff_fourier'];
path='C:\Users\manon\Desktop\projet CSE I\figures\cas_ecole';


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