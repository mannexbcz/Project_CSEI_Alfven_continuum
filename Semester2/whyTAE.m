% This scripts compares the values of M,N, GradPsi2, R^2, B^2 for different
% triangularites (Figure 22)
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%%
save=1; %Set 1 to save the figure, 0 otherwise
path='C:\Users\manon\Desktop\projet CSE II\figures\fig_whyTAE\';

%%

folder = uigetdir; %gets directory
fileList = dir(fullfile(folder, '*.h5'));
qval = 1.5;

Ms = [];
Ns = [];
gradpsis=[];
B2s = [];
R2s = [];
deltas = [];
CoeffsMs = [];
CoeffsNs = [];

for k = 1:3:length(fileList)
    
   
    baseFileName = fileList(k).name;
    fullFileName = fullfile(folder, baseFileName);
    params = read_chease(fullFileName);
    
    a = params.a; R0 = params.R0; B0 = params.B0;
    q = params.q; B = params.B; R = params.R; F = params.F;
    gradPsi2 = params.gradPsi2;
    r = params.r;
    nr = length(r);
    thetastar = (params.thetastar)';
    deltas = [deltas,params.delta(end)];
    
    
    M_all = gradPsi2.*(B0./B).^2.*(R./R0).^2;
    N_all = gradPsi2./(B.^2.*R.^2);
    
    idx = index_q_r(q,qval);
    
    Ms = [Ms, M_all(idx,:)'];
    Ns = [Ns, N_all(idx,:)'];
    gradpsis = [gradpsis, gradPsi2(idx,:)'];
    B2s = [B2s, B(idx,:)'.^2];
    R2s = [R2s, R(idx,:)'.^2];
    
    %Fourier Coefficients
    
    band = 11;
    m_min_coeff=0;
    m_max_coeff=(band+1)/2;
    
   % Fourier coefficients(fun,thetastar,mmin,mmax)
    [coeffs_M] = get_fourier_coeff_chease(M_all(idx,:),thetastar,m_min_coeff,m_max_coeff);
    [coeffs_N] = get_fourier_coeff_chease(N_all(idx,:),thetastar,m_min_coeff,m_max_coeff);
    
    CoeffsMs = [CoeffsMs, coeffs_M'];
    CoeffsNs = [CoeffsNs, coeffs_N'];

    
end

%% Figure M

figure
plot(thetastar, Ms)
legend(strcat('$\delta_a= $',num2str(deltas(1),'%.2f')),strcat('$\delta_a= $',num2str(deltas(2),'%.2f')),strcat('$\delta_a= $',num2str(deltas(3),'%.2f')),...
    strcat('$\delta_a= $',num2str(deltas(4),'%.2f')),strcat('$\delta_a= $',num2str(deltas(5),'%.2f')))
grid on
%title('$q(r)=1.5$')
xlabel('$\theta^*$')
ylabel('$M$')

namefig = 'M_diff_triang';
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

%% Figure N

figure
plot(thetastar, Ns)
legend(strcat('$\delta_a= $',num2str(deltas(1),'%.2f')),strcat('$\delta_a= $',num2str(deltas(2),'%.2f')),strcat('$\delta_a= $',num2str(deltas(3),'%.2f')),...
    strcat('$\delta_a= $',num2str(deltas(4),'%.2f')),strcat('$\delta_a= $',num2str(deltas(5),'%.2f')))
grid on
%title('$q(r)=1.5$')
xlabel('$\theta^*$')
ylabel('$N$')

namefig = 'N_diff_triang';
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
%% Figure GradPsi2

figure
plot(thetastar, gradpsis)
legend(strcat('$\delta_a= $',num2str(deltas(1),'%.2f')),strcat('$\delta_a= $',num2str(deltas(2),'%.2f')),strcat('$\delta_a= $',num2str(deltas(3),'%.2f')),...
    strcat('$\delta_a= $',num2str(deltas(4),'%.2f')),strcat('$\delta_a= $',num2str(deltas(5),'%.2f')))
grid on
%title('$q(r)=1.5$')
xlabel('$\theta^*$')
ylabel('$\nabla \Psi ^2$')

namefig = 'GradPsi2_diff_triang';
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
%% Figure B^2

figure
plot(thetastar, B2s)
legend(strcat('$\delta_a= $',num2str(deltas(1),'%.2f')),strcat('$\delta_a= $',num2str(deltas(2),'%.2f')),strcat('$\delta_a= $',num2str(deltas(3),'%.2f')),...
    strcat('$\delta_a= $',num2str(deltas(4),'%.2f')),strcat('$\delta_a= $',num2str(deltas(5),'%.2f')))
grid on
%title('$q(r)=1.5$')
xlabel('$\theta^*$')
ylabel('$B^2$')

namefig = 'B2_diff_triang';
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
%% Figure R^2

figure
plot(thetastar, R2s)
legend(strcat('$\delta_a= $',num2str(deltas(1),'%.2f')),strcat('$\delta_a= $',num2str(deltas(2),'%.2f')),strcat('$\delta_a= $',num2str(deltas(3),'%.2f')),...
    strcat('$\delta_a= $',num2str(deltas(4),'%.2f')),strcat('$\delta_a= $',num2str(deltas(5),'%.2f')))
grid on
%title('$q(r)=1.5$')
xlabel('$\theta^*$')
ylabel('$R^2$')

namefig = 'R2_diff_triang';
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

%% Figure coeffs
x=[1,2,3,4,5,6,7];
figure
plot(x,CoeffsMs)
legend(strcat('$\delta_a= $',num2str(deltas(1),'%.2f')),strcat('$\delta_a= $',num2str(deltas(2),'%.2f')),strcat('$\delta_a= $',num2str(deltas(3),'%.2f')),...
    strcat('$\delta_a= $',num2str(deltas(4),'%.2f')),strcat('$\delta_a= $',num2str(deltas(5),'%.2f')))
grid on
%title('$q(r)=1.5$')
xlabel('Number of the coefficient')
ylabel('Fourier coefficient of M')

namefig = 'Mcoeffs_diff_triang';
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

figure
plot(x,CoeffsNs)
legend(strcat('$\delta_a= $',num2str(deltas(1),'%.2f')),strcat('$\delta_a= $',num2str(deltas(2),'%.2f')),strcat('$\delta_a= $',num2str(deltas(3),'%.2f')),...
    strcat('$\delta_a= $',num2str(deltas(4),'%.2f')),strcat('$\delta_a= $',num2str(deltas(5),'%.2f')))
grid on
%title('$q(r)=1.5$')
xlabel('Number of the coefficient')
ylabel('Fourier coefficient of N')

namefig = 'Ncoeffs_diff_triang';
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end