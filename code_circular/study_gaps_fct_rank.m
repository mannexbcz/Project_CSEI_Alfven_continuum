% -------------------------------------------------------------------------
% This script explores the impact on the convergence of the rank of the 
% matrices M and N considered. A circular, concentric and large aspect
% ratio equilibrium is considered.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters

% rank of the matrices M and N
sizes=[2:15];  

a=0.1;      % minor radius
n=1;        % toroidal mode number
nr = 1000;  % number of magnetic surfaces considered
R0=1;       % major radius

q = @(r) 1+2*(r/a).^2; % safety factor


%% Scan over the rank

gaps = zeros(1,length(sizes));

for s=1:length(sizes)

size=sizes(s);
    
w = eigenmodes_circular(size,a,n,R0,nr,q);

% Detects peaks of the first and second eigenmodes, associated with TAE gaps.
[pks1,locs1] = findpeaks(w(1,:),r./a,'MinPeakDistance',0.1);
[pks2,locs2] = findpeaks(-w(2,:),r./a,'MinPeakDistance',0.1);

% TAE gap
gaps(s)=-pks2(1)-pks1(1);

end

%% Reference

size=100;

w = eigenmodes_circular(size,a,n,R0,nr,q);

[pks1,locs1] = findpeaks(w(1,:),r./a,'MinPeakDistance',0.1);
[pks2,locs2] = findpeaks(-w(2,:),r./a,'MinPeakDistance',0.1);

gaps_ref=-pks2(1)-pks1(1);

%% Figure

Const1=polyfit(log(sizes(1:8)),log(abs(gaps(1:8)-gaps_ref)./gaps_ref),1);
m = Const1(1);
k= Const1(2);
JBL=10.^(m*sizes(1:8) + k);
JBL = sizes(1:8).^m.*exp(k);

save=1; %Set 1 to save the figure, 0 otherwise
% name = ['Rank of A,B =', num2str(size)] ;
namefig=['convergence_gap_rank'];
path='C:\Users\manon\Desktop\projet CSE I\figures\circulaire\';


figure
loglog(sizes,abs(gaps-gaps_ref)./gaps_ref,'k +','HandleVisibility','off')
hold on
% loglog(sizes(1:8),JBL)
grid on
xlabel('Rank')
ylabel('Relative error on the Gap Size')
% legend('$y=-19.1x + b$')
% title(name)

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end