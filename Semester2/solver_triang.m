% -------------------------------------------------------------------------
% This code solves the matrix form of the Alfvén continuum equation and 
% provides the Alfvén spectrum.
% 4 Different Options:
%   - Geometric : Solution for simplified equilibria parametrized using 
%   Miller's formulas
%   - CHEASEGeom : solution via geometrical parametrization and 
%   coefficients obtained from CHEASE
%   - CHEASE : Solution without geometrical parametrization using 
%   equilibrium coefficients directly obtained from CHEASE
%   - Quadratic : Solution with quadratic approximation of the q-profile
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
n=1;         % toroidal mode number
npoints=512; % number of mesh points used for the integration

Option = input('Method (Geometric/CHEASE/CHEASEGeom/Quadratic): ','s');

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
    q = @(r) 1+2*(r/a).^2; 
    kfct = @(r) kprime.*(r-a) + k_a;
    deltafct = @(r)  -deltaprime.*a + deltaprime.*r;
    dfct = @(r) dprime.*r;
    
    % Eigenmodes
    w = eigenmodes_triang(size,band,kfct,k_a,kprime,deltafct,deltaprime,dfct,dprime,q,R0,B0,a,n,nr,npoints);

elseif strcmp(Option,'CHEASE')
    [filename, pathname] = uigetfile( {'*.*','All Files (*.*)'},'Pick the CHEASE equilibria file');
    if(filename==0)
        warning('No file selected.')
    end
    path = append(pathname,filename);
    params = read_chease(path);
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
    r = params.r;
    nr = length(r);
    a = params.a;
    w = eigenmodes_chease_geom(size,band,params,n,nr,npoints);

elseif strcmp(Option,'Quadratic')
    [filename, pathname] = uigetfile( {'*.*','All Files (*.*)'},'Pick the CHEASE equilibria file');
    if(filename==0)
        warning('No file selected.')
    end    
    params = read_chease(filename);
    r = params.r;
    nr = length(r);
    a = params.a; 
    
    w = eigenmodes_quadratic_q(size,band,params,n,nr,npoints);
    
else
    error('Invalid Method')
end

%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
namefig = 'triang_pos';
path='C:\Users\manon\Desktop\projet CSE I\figures\triang\';

figure
for i=1:size-1
plot(r./a,w(i,:),'b', 'HandleVisibility','off')
hold on
end
plot(r./a,w(end,:),'b')
grid on
xlabel('$r/a$')
ylabel('$\tilde{\omega}$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
 
