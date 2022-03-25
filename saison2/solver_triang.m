% -------------------------------------------------------------------------
% This code solves the matrix form of the Alfvén continuum equation in the 
% case of an elongated and shifted equilibrium, with triangularity.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters
clear all

size=15;     % rank of the matrices M and N
band=11;     % bands of the matrices M and N
n=1;        % toroidal mode number
npoints=512;% number of mesh points used for the integration

Option = input('Method (Geometric/CHEASE/CHEASEGeom): ','s');

if strcmp(Option,'Geometric')
    a=0.1;      % minor radius
    R0=1;       % major radius
    B0=1;       % magnetic field
    F = R0*B0;  % toroidal flux current
    % Elongation
    k_a=1.1;
    kprime=0;
    %Shift
    deltaprime= -0.1;
    %Triangularity
    dprime = -0.5;
    nr = 301;  % number of magnetic surfaces considered

    % Functions
    kfct = @(r) kprime.*(r-a) + k_a;
    deltafct = @(r)  -deltaprime.*a + deltaprime.*r;
    dfct = @(r) dprime.*r;
    
    % Eigenmodes
    % r=transpose(linspace(0,a,nr));
    w = eigenmodes_triang(size,band,kfct,k_a,kprime,deltafct,deltaprime,dfct,dprime,q,R0,B0,a,n,nr,npoints,Option);

elseif strcmp(Option,'CHEASE')
    [filename, pathname] = uigetfile( {'*.*','All Files (*.*)'},'Pick the CHEASE equilibria file');
    if(filename==0)
        warning('No file selected.')
    end
    params = read_chease(filename);
    % [a,R0,B0,k_a,kprime,deltaprime,dprime,q]
    a = params.a; R0 = params.R0; B0 = params.B0;
    q = params.q; B = params.B; R = params.R; F = params.F;
    gradPsi2 = params.gradPsi2;
    r = params.r;
    nr = length(r);
    thetastar = (params.thetastar)';
 
    w = eigenmodes_chease(size,band,gradPsi2,B,B0,R,R0,F,q,thetastar,a,n,nr);

elseif strcmp(Option,'CHEASEGeom')
    [filename, pathname] = uigetfile( {'*.*','All Files (*.*)'},'Pick the CHEASE equilibria file');
    if(filename==0)
        warning('No file selected.')
    end    
    params = read_chease(filename);
    % [a,R0,B0,k_a,kprime,deltaprime,dprime,q]
    r = params.r;
    nr = length(r);
    a = params.a;
    w = eigenmodes_chease_geom(size,band,params,n,nr,npoints);
    
else
    error('Invalid Method')
end

%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
namefig = 'triang_pos';
path='C:\Users\manon\Desktop\projet CSE I\figures\triang\';

figure
for i=1:size-1
plot(r./a,w(i,:),'k', 'HandleVisibility','off')
hold on
end
plot(r./a,w(end,:),'k')
grid on
xlabel('$r/a$')
ylabel('$\tilde{\omega}$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
 
