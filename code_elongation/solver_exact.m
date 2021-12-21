% -------------------------------------------------------------------------
% This code solves the matrix form of the Alfvén continuum equation in the 
% case of an elongated and shifted equilibrium
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters

size=3;     % rank of the matrices M and N
band=5;     % bands of the matrices M and N
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
B0=1;       % magnetic field
F = R0*B0;  % toroidal flux current
nr = 100;  % number of magnetic surfaces considered
npoints=100;% number of mesh points used for the integration
q = @(r) 1+2*(r/a).^2; % safety factor

% Elongation
k_a=1;
kprime=0;
kfct = @(r) kprime.*(r-a) + k_a;

%Shift
deltaprime= 0;
deltafct = @(r)  -deltaprime*a + deltaprime *r;

%% Eigenmodes
r=linspace(0,a,nr);
w = eigenmodes_exact(size,band,kfct,k_a,kprime,deltafct,deltaprime,q,R0,B0,a,n,nr,npoints);


%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
name = ['Rank of A,B =', num2str(size)] ;
namefig=['frequences_rank=',num2str(size)];
path='C:\Users\manon\Desktop\projet CSE I\figures\';

figure
for i=1:size-1
plot(r./a,w(i,:),'r-','HandleVisibility', 'off')
hold on
end
plot(r./a,w(end,:),'r-')
grid on
xlabel('$r/a$')
ylabel('$\tilde{\omega}$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
 
