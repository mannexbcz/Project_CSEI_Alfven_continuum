% -------------------------------------------------------------------------
% This script studies the impact of the triangularity on the NAE and TAE gaps 
% An elongated and shifted equilibrium, with triangularity, is used.
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
band=7;     % bands of the matrices M and N
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
B0=1;       % magnetic field
F = R0*B0;  % toroidal flux current
nr = 100;  % number of magnetic surfaces considered
npoints=100;% number of mesh points used for the integration
q = @(r) 1+2*(r/a).^2; % safety factor

% Elongation
k_a=1.1;
kprime=0;
kfct = @(r) kprime.*(r-a) + k_a;

%Shift
deltaprime= -0.1;
deltafct = @(r)  -deltaprime*a + deltaprime *r;

%Triangularity
dprimes = linspace(-5,5,30);

%%
TAEgaps=zeros(1,length(dprimes));
EAEgaps=zeros(1,length(dprimes));
NAEgaps=zeros(1,length(dprimes));
r=linspace(0,a,nr);

for s = 1:length(dprimes)
    dprime = dprimes(s);
    dfct = @(r) dprime.*r;
    w = eigenmodes_triang(size,band,kfct,k_a,kprime,deltafct,deltaprime,dfct,dprime,q,R0,B0,a,n,nr,npoints);

    %TAE gaps
    [pks1,locs1] = findpeaks(real(w(1,:)),'MinPeakDistance',10);
    [pks2,locs2] = findpeaks(-real(w(2,:)),'MinPeakDistance',10);

    TAEgaps(s)=w(2,locs2(1))-w(1,locs1(1));
    
    % EAE gaps
    [pks3,locs3] = findpeaks(real(w(2,:)));
    [pks4,locs4] = findpeaks(-real(w(3,:))); 
    
    EAEgaps(s)=-max(pks4)-min(pks3); %w(3,locs4(1))-w(2,locs3(1));

    % NAE gaps
    [pks5,locs5] = findpeaks(real(w(3,:)));
    [pks6,locs6] = findpeaks(-real(w(4,:))); 

    NAEgaps(s)=-max(pks6)-min(pks4); 
end
%% Figures

figure
plot(dprimes.*a,TAEgaps,'k +')
grid on
xlabel('$\delta_a$')
ylabel('TAE Gap')

figure
plot(dprimes.*a,EAEgaps,'k +')
grid on
xlabel('$\delta_a$')
ylabel('EAE Gap')

figure
plot(dprimes.*a,NAEgaps,'k +')
grid on
xlabel('$\delta_a$')
ylabel('NAE Gap')