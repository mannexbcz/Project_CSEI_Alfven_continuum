% -------------------------------------------------------------------------
% This code solves the matrix form of the Alfvén continuum equation in the 
% textbook case (eq.(78))
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
nr = 1000;          % number of magnetic surfaces considered
npoints = 1000;     % number of mesh points for the numerical integration
q = @(r) 1+2*(r/a).^2; % safety factor

%% Equilibrium coefficients
M = @(epsilon,epsilon2,epsilon3,theta) 1+4.*epsilon.*cos(theta+epsilon2.*sin(theta));
N = @(epsilon,epsilon2,epsilon3,theta) 1+epsilon3.*cos(theta+epsilon2.*sin(theta)).^2;

epsilon = 0.1;
epsilon2 = 0.2;
epsilon3 = 0.3;

%% Eigenmodes
r=linspace(0,a,nr);
w = eigenmodes_textbookcase(M,N,size,band,epsilon,epsilon2,epsilon3,a,n,F,nr,q,npoints);

%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['frequences_rank=',num2str(size)];
path='C:\Users\manon\Desktop\projet CSE I\figures\';

figure
for i=1:size
plot(r,w(i,:))
hold on
end
grid on
xlabel('$r$')
ylabel('$\omega$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end