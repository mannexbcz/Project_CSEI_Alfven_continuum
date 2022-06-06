% -------------------------------------------------------------------------
% This code solves the matrix form of the Alfvén continuum equation in the 
% case of a circular and shifted equilibrium. The matrix form is given by 
% eq. (19) of the report, with matrices M,N given by eq. (39),(40).
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1.5);
format long;

%% Parameters

size=3;     % rank of the matrices M and N
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
B0=1;       % magnetic field
delta_p=-0.0;% derivative of the shift
F = R0*B0;  % toroidal flux current
nr = 201;  % number of magnetic surfaces considered
q = @(r) 1+2*(r/a).^2; % safety factor

%% Alfvén eigenmodes

r=linspace(0,a,nr); % magnetic surfaces considered
w = eigenmodes_shifted(size,delta_p,a,n,R0,B0,nr,q); % Alfvén eigenmodes

%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
name = ['Rank of A,B =', num2str(size)] ;
namefig=['solver_shifted'];
path='C:\Users\manon\Desktop\projet CSE I\figures\';

figure
plot(r./a,w(1,:),'k')
hold on
plot(r./a,w(2,:),'k','HandleVisibility','off')
hold on
plot(r./a,w(3,:),'k','HandleVisibility','off')
grid on
xlabel('$r/a$')
ylabel('$\omega$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
