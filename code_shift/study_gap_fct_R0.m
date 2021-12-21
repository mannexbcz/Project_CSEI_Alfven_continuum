% -------------------------------------------------------------------------
% This script explores the impact of the inverse aspect ratio on the size
% of the first TAE gap, in the case of a circular but shifted equilibriuum.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Paramètres

sizes=[2:15];     % rank of the matrices M and N
a=0.1;      % minor radius
n=1;        % toroidal mode number
B0=1;       % magnetic field
delta_p=-0.1;% derivative of the shift
F = R0*B0;  % toroidal flux current
nr = 1000;  % number of magnetic surfaces considered
q = @(r) 1+2*(r/a).^2; % safety factor

R0s=logspace(-1,1,40); % major radius

%% Construction des matrices & résolution du système

gaps=zeros(1,length(R0s));

r=linspace(0,a,1000);

for s=1:length(R0s)

R0=R0s(s);
    
w = eigenmodes_shifted(size,delta_p,a,n,R0,B0,nr,q); % Alfvén eigenmodes

[pks1,locs1] = findpeaks(w(1,:),'MinPeakDistance',100);
[pks2,locs2] = findpeaks(-w(2,:),'MinPeakDistance',100);

gaps(s)=w(2,locs1(1))-w(1,locs1(1));

end

%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['gap_fct_R0_shifted'];
path='C:\Users\manon\Desktop\projet CSE I\figures\shifted\';

figure
semilogy(a./R0s,gaps,'k +')
grid on
xlabel('$a/R_0$')
ylabel('Gap Size')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end