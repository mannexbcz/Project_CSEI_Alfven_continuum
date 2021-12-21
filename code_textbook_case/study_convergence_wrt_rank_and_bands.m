% -------------------------------------------------------------------------
% This script explores the impact on the convergence of the rank and the
% number of bands of the matrices M and N considered, in the textbook case
% studied (eq.(78)).
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters

a=0.1;              % minor radius
n=1;                % toroidal mode number
R0=1;               % major radius
B0=1;               % magnetic field
F = R0*B0;          % toroidal flux current
nr = 100;          % number of magnetic surfaces considered
npoints = 1000;     % number of mesh points for the numerical integration
q = @(r) 1+2*(r/a).^2; % safety factor

%% Equilibrium coefficients

M = @(epsilon,epsilon2,epsilon3,theta) 1+4.*epsilon.*cos(theta+epsilon2.*sin(theta));
N = @(epsilon,epsilon2,epsilon3,theta) 1+epsilon3.*cos(theta+epsilon2.*sin(theta)).^2;

epsilon = 0.1;
epsilon2 = 0.2;
epsilon3 = 0.3;

%% Convergence with respect to the number of bands of the matrices

size=10;             % rank of the matrices M and N
bands=[3,5,7,9,11,13,15,17,19];

gaps=zeros(1,length(bands));

for s=1:length(bands)

band=bands(s);
w = eigenmodes_textbookcase(M,N,size,band,epsilon,epsilon2,epsilon3,a,n,F,nr,q,npoints);

[pks1,locs1] = findpeaks(real(w(1,:)),'MinPeakDistance',10);
[pks2,locs2] = findpeaks(-real(w(2,:)),'MinPeakDistance',10);

gaps(s)=w(2,locs1(1))-w(1,locs1(1));

end

%% Figure

save=1; %Set 1 to save the figure, 0 otherwise
namefig=['cas_ecole_gap_fct_band'];
path='C:\Users\manon\Desktop\projet CSE I\figures\cas_ecole\';
gaps=real(gaps);
diff=gaps-gaps(end);

figure
semilogy(bands,abs(diff)./gaps,'k +')
grid on
xlabel('Bands')
ylabel('Relative error on the Gap Size')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

%% Convergence with respect to the rank of the matrices

band = 3;
sizes=[2:15]; % rank of the matrices
gaps=zeros(1,length(sizes));

for s=1:length(sizes)

size=sizes(s);
w = eigenmodes_textbookcase(M,N,size,band,epsilon,epsilon2,epsilon3,a,n,F,nr,q,npoints);

[pks1,locs1] = findpeaks(real(w(1,:)),'MinPeakDistance',10);
[pks2,locs2] = findpeaks(-real(w(2,:)),'MinPeakDistance',10);

gaps(s)=w(2,locs1(1))-w(1,locs1(1));

end

%% Reference

size=100;
w = eigenmodes_textbookcase(M,N,size,band,epsilon,epsilon2,epsilon3,a,n,F,nr,q,npoints);

[pks1,locs1] = findpeaks(real(w(1,:)),'MinPeakDistance',10);
[pks2,locs2] = findpeaks(-real(w(2,:)),'MinPeakDistance',10);

gaps_ref=w(2,locs1(1))-w(1,locs1(1));

%% Figure

save=1; %Set 1 to save the figure, 0 otherwise
namefig=['cas_ecole_gap_fct_rank'];
path='C:\Users\manon\Desktop\projet CSE I\figures\cas_ecole\';

figure
semilogy(sizes,abs((gaps-gaps_ref)./gaps_ref),'k +')
grid on
xlabel('Rank')
ylabel('Relative error on the Gap Size')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end