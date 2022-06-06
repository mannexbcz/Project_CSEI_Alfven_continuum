% -------------------------------------------------------------------------
% This script studies the impact on the convergence of the rank and the
% number of bands of the matrices M and N considered, for the solution
% using the equilibrium coefficients from CHEASE
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters
n=1;        % toroidal mode number
npoints=512;% number of mesh points used for the integration

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

%% Convergence with respect to the number of bands of the matrices
size=15;     % rank of the matrices M and N
bands=[3,5,7,9,11,13,15,17,19,21,23,25,27,29];  % number of bands of the matrices M and N

gaps=zeros(1,length(bands));
time_bands = zeros(1,length(bands));

r=linspace(0,a,nr);

for s=1:length(bands)

band=bands(s);

w = eigenmodes_chease(size,band,gradPsi2,B,B0,R,R0,F,q,thetastar,a,n,nr);

[TAE,reltativeTAE] = gap_at_q(1.5,q,w(1,:),w(2,:));

gaps(s)=TAE;

end

%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['ordre_1_error_gap_fct_nb_bandes'];
path='C:\Users\manon\Desktop\projet CSE I\figures\elongation\';
gaps=real(gaps);
diff=gaps-gaps(end);

figure
semilogy(bands,abs(diff)/gaps(end)./gaps,'k +')
grid on
xlabel('Bands')
ylabel('Error on the Gap Size')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

%% Convergence with respect to the rank of the matrices

sizes=[2:40];  % rank of the matrices M and N
band = 11; % number of bands of the matrices M and N

gaps=zeros(1,length(sizes));
time_size = zeros(1,length(sizes));

for s=1:length(sizes)

size=sizes(s);
    
w = eigenmodes_chease(size,band,gradPsi2,B,B0,R,R0,F,q,thetastar,a,n,nr);

[TAE,reltativeTAE] = gap_at_q(1.5,q,w(1,:),w(2,:));

gaps(s)=TAE;
end

%% Reference

size=100;
band = 11;

w = eigenmodes_chease(size,band,gradPsi2,B,B0,R,R0,F,q,thetastar,a,n,nr);

[TAE,reltativeTAE] = gap_at_q(1.5,q,w(1,:),w(2,:));

gaps_ref=TAE;

%% Figure

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['ordre_1_error_gap_fct_rang'];
path='C:\Users\manon\Desktop\projet CSE I\figures\elongation\';

figure
semilogy(sizes,abs((gaps-gaps_ref)./gaps_ref),'k +')
grid on
xlabel('Rank')
ylabel('Error on the Gap Size')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end