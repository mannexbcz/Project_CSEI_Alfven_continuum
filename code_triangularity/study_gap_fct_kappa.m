% -------------------------------------------------------------------------
% This script studies the impact on the EAE gap size of the value of k_a. 
% An elongated and shifted equilibrium is used, in the first order
% approximation.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters
size=5;     % rank of the matrices M and N
band=5;     % bands of the matrices M and N
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
B0=1;       % magnetic field
F = R0*B0;  % toroidal flux current
nr = 100;  % number of magnetic surfaces considered
npoints=100;% number of mesh points used for the integration
q = @(r) 1+2*(r/a).^2; % safety factor

%Shift
deltaprime= -0.1;
deltafct = @(r)  -deltaprime*a + deltaprime *r;
%Triangularity
dprime = 0;
dfct = @(r) dprime.*r;

kappas=linspace(1,1.7,20);

%%
gaps=zeros(1,length(kappas));
r=linspace(0,a,nr);

for s = 1:length(kappas)
    k_a=kappas(s);
    kprime=0; 
    kfct = @(r) kprime.*(r-a) + k_a;
    
    w = eigenmodes_triang(size,band,kfct,k_a,kprime,deltafct,deltaprime,dfct,dprime,q,R0,B0,a,n,nr,npoints);

% [pks1,locs1] = findpeaks(real(w(2,50:end)));
% [pks2,locs2] = findpeaks(-real(w(3,50:end))); 
diff = abs(real(w(3,50:75))-real(w(2,50:75)));

gaps(s)=min(diff);

end

%% Figure

save=1; %Set 1 to save the figure, 0 otherwise
namefig=['gap_fct_kappa'];
path='C:\Users\manon\Desktop\projet CSE I\figures\elongation\';


figure
plot(kappas,gaps,'k +')
grid on
xlabel('$\kappa_a$')
ylabel('Gap Size')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
