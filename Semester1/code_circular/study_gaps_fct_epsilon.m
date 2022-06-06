% -------------------------------------------------------------------------
% This script explores the impact of the inverse aspect ratio on the size
% of the first TAE gap.
% It solves the Alfvén continuum equation in the case of a circular, 
% concentric, and large aspect equilibriuum, and for different values of 
% the major radius R0, that is for different values of the inverse aspect 
% ratio epsilon. The first TAE gap is then deduced from the obtained 
% eigenmodes.
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
nr = 1000;  % number of magnetic surfaces considered
q = @(r) 1+2*(r/a).^2; % safety factor

% Values of the major radius R0 considered
R0s = logspace(-1,1,40);

%% Scan over the values of R0
% For each value of R0, the Alfvén equation is solved

gaps=zeros(1,length(R0s));
r=linspace(0,a,nr);

for s=1:length(R0s) 

R0=R0s(s); 
  
w = eigenmodes_circular(size,a,n,R0,nr,q);

% Detects peaks of the first and second eigenmodes, associated with TAE gaps.
[pks1,locs1] = findpeaks(w(1,:),'MinPeakDistance',100);
[pks2,locs2] = findpeaks(-w(2,:),'MinPeakDistance',100);

% First TAE gap
gaps(s)=w(2,locs1(1))-w(1,locs1(1));

end

%% Figure
% Plots the size of the first TAE gap as a function of R0

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['gap_fct_R0'];
path='C:\Users\manon\Desktop\projet CSE I\figures\circulaire\';


figure
semilogy(a./R0s,gaps,'k +')
grid on
xlabel('$a/R_0$')
ylabel('Gap Size')
% title(name)

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
