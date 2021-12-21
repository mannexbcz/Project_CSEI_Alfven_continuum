% -------------------------------------------------------------------------
% This code solves the matrix form of the Alfvén continuum equation in the 
% case of a circular, concentric, and large aspect ratio equilibrium. The
% matrix form is given by eq. (19) of the report, with matrices M,N given
% by eq. (28),(29).
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters

size=2;     % rank of the matrices M and N
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
B0=1;       % magnetic field
nr = 1000;  % number of magnetic surfaces considered
q = @(r) 1+2*(r/a).^2; % safety factor

%% Alfvén eigenmodes

r=linspace(0,a,nr); % magnetic surfaces considered
w = eigenmodes_circular(size,a,n,R0,B0,nr,q); % Alfvén eigenmodes

%% Figure
% This figure plots the obtained Alfvén contiuum

save=0; %Set 1 to save the figure, 0 otherwise
name = ['Rank of A,B =', num2str(size)] ;
namefig=['frequences_rank=',num2str(size)];
path='C:\Users\manon\Desktop\projet CSE I\figures\circulaire\';


figure
for i=1:size
plot(r./a,w(i,:))
hold on
end

grid on
xlabel('$r/a$')
ylabel('$\tilde{\omega}$')
%title(name)

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end