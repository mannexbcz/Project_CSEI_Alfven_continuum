% -------------------------------------------------------------------------
% This script studies the impact on the gap size of the value of delta prime. 
% A circular and shifted equilibrium is used.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1.5);
format long;

%% Parameters

size=2;     % rank of the matrices M and N
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
B0=1;       % magnetic field
F = R0*B0;  % toroidal flux current
nr = 1000;  % number of magnetic surfaces considered
q = @(r) 1+2*(r/a).^2; % safety factor

delta_ps=linspace(0,0.8,30); % derivative of the shift

%% Scan over delta prime

gaps=zeros(1,length(delta_ps));

for s=1:length(delta_ps)

delta_p=delta_ps(s);

w = eigenmodes_shifted(size,delta_p,a,n,R0,B0,nr,q); % Alfvén eigenmodes

[pks1,locs1] = findpeaks(w(1,:),'MinPeakDistance',100);
[pks2,locs2] = findpeaks(-w(2,:),'MinPeakDistance',100);

gaps(s)=w(2,locs1(1))-w(1,locs1(1));

end

%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['gap_fct_deltap'];
path='C:\Users\manon\Desktop\projet CSE I\figures\shifted\';


figure
plot(delta_ps,gaps,'k +')
grid on
xlabel('$\Delta^\prime  $')
ylabel('Gap Size')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
