% -------------------------------------------------------------------------
% This function solves the matrix form of the Alfvén continuum equation and 
% provides the Alfvén spectrum.
% -------------------------------------------------------------------------
%%
function solver()

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
end
