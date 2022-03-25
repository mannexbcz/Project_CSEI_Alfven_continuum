% -------------------------------------------------------------------------
% This code solves the matrix form of the Alfvén continuum equation in the 
% case of an elongated and shifted equilibrium, with triangularity.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters

size=10;     % rank of the matrices M and N
band=7;     % bands of the matrices M and N
a=0.24;      % minor radius
n=1;        % toroidal mode number
R0=0.88;       % major radius
B0=1;       % magnetic field
F = R0*B0;  % toroidal flux current
nr = 201;  % number of magnetic surfaces considered
npoints=33;% number of mesh points used for the integration
q = @(r) 1+2*(r/a).^2; % safety factor

% Elongation
k_a=1.4;
kprime=0;
kfct = @(r) kprime.*(r-a) + k_a;

%Shift
deltaprime= -0.1;
deltafct = @(r)  -deltaprime*a + deltaprime *r;

%Triangularity
dprime = 0.8/a;
dfct = @(r) dprime.*r;

%% Eigenmodes
tic
r=linspace(0,a,nr);
w = eigenmodes_triang(size,band,kfct,k_a,kprime,deltafct,deltaprime,dfct,dprime,q,R0,B0,a,n,nr,npoints);
toc

%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
namefig = 'triang_pos';
path='C:\Users\manon\Desktop\projet CSE I\figures\triang\';

figure
for i=1:size-1
plot(r./a,w(i,:),'k', 'HandleVisibility','off')
hold on
end
plot(r./a,w(end,:),'k')
grid on
xlabel('$r/a$')
ylabel('$\tilde{\omega}$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
 
