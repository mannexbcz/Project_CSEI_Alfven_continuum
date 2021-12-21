% -------------------------------------------------------------------------
% This script explores the impact on the convergence of the rank of the 
% matrices M and N considered. A circular and shifted equilibrium is used.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters

sizes=[2:15];     % rank of the matrices M and N
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
B0=1;       % magnetic field
delta_p=-0.1;% derivative of the shift
F = R0*B0;  % toroidal flux current
nr = 1000;  % number of magnetic surfaces considered
q = @(r) 1+2*(r/a).^2; % safety factor

%% Scan over the rank

gaps=zeros(1,length(sizes));

for s=1:length(sizes)

size=sizes(s);
    
w = eigenmodes_shifted(size,delta_p,a,n,R0,B0,nr,q); % Alfvén eigenmodes

% Detects peaks of the first and second eigenmodes, associated with TAE gaps.
[pks1,locs1] = findpeaks(w(1,50:end),r(50:end)./a,'MinPeakDistance',0.1);
[pks2,locs2] = findpeaks(-w(2,50:end),r(50:end)./a,'MinPeakDistance',0.1);

% TAE gap
gaps(s)=-pks2(1)-pks1(1);

end

%% Reference

size=100;
w = eigenmodes_shifted(size,delta_p,a,n,R0,B0,nr,q); % Alfvén eigenmodes

[pks1,locs1] = findpeaks(w(1,50:end),r(50:end)./a,'MinPeakDistance',0.1);
[pks2,locs2] = findpeaks(-w(2,50:end),r(50:end)./a,'MinPeakDistance',0.1);

gaps_ref=-pks2(1)-pks1(1);

%% figure 2

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['shifted_error_conv_gap_fct_rank'];
path='C:\Users\manon\Desktop\projet CSE I\figures\shifted\';

figure
semilogy(sizes,abs((gaps-gaps_ref)./gaps_ref),'k +')
grid on
xlabel('Rank')
ylabel('Relative error on the Gap Size')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end