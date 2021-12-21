% -------------------------------------------------------------------------
% This script explores the impact on the convergence of the number of 
% magnetic surfaces considered. A circular, concentric and large aspect
% ratio is considered.
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
q = @(r) 1+2*(r/a).^2; % safety factor

rmaxs=logspace(1,4,20); % number of magnetic surfaces considered
%% Scan over the number of magnetic surfaces considered

gaps=zeros(1,length(rmaxs));

for s=1:length(rmaxs)
    
nr=rmaxs(s);

w = eigenmodes_circular(size,a,n,R0,round(nr),q);

% Detects peaks of the first and second eigenmodes, associated with TAE gaps.
[pks1,locs1] = findpeaks(w(1,:),'MinPeakProminence',0.05);
[pks2,locs2] = findpeaks(-w(2,:),'MinPeakProminence',0.05);

% First TAE gap
gaps(s)=w(2,locs2(1))-w(1,locs1(1));

end
%% Reference

nr = 5e4;    
w = eigenmodes_circular(size,a,n,R0,round(nr),q);


% Detects peaks of the first and second eigenmodes, associated with TAE gaps.
[pks1,locs1] = findpeaks(w(1,:),'MinPeakProminence',0.1);
[pks2,locs2] = findpeaks(-w(2,:),'MinPeakProminence',0.1);

% Reference for the first TAE gap
gaps_ref=w(2,locs2(1))-w(1,locs1(1));


%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['error_gap_fct_N'];
path='C:\Users\manon\Desktop\projet CSE I\figures\circulaire\';

Const1=polyfit(log(rmaxs),log(abs(gaps-gaps_ref)./gaps_ref),1);
m = Const1(1);
k= Const1(2);
JBL=10.^(m*rmaxs + k);
JBL = rmaxs.^m.*exp(k);

figure
loglog(rmaxs,abs((gaps-gaps_ref)./gaps_ref),'k +','HandleVisibility','off')
hold on
loglog(rmaxs,JBL)
grid on
xlabel('$N$')
ylabel('Relative error on the Gap Size')
legend('$y=-1.9x+b$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end