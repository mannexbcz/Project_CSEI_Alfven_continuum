% This script is used to plot the values of the gap for different
% triangularities. It reads the outputs of CHEASE obtained for various 
% triangularities, computes the Alfvén spectrum and corresponding gaps and
% plots the results.

%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%%
size=20;     % rank of the matrices M and N
band=15;     % bands of the matrices M and N
npoints=512;% number of mesh points used for the integration

folder = uigetdir; %gets directory

fileList = dir(fullfile(folder, '*.h5'));
%%
n=2;        % toroidal mode number

TAEgaps=[];
EAEgaps=[];
NAEgaps=[];
deltas =[];
deltasTAE=[];
deltasEAE=[];
deltasNAE=[];

switch n
    case 1
        val1 = 1.5; val2=2; val3=2.5;
        c= 'b-+';
    case 2
        val1 = 5/4; val2=3/2; val3=7/4;
        c= 'r-+';
    case 3
        val1 = 7/6; val2=4/3; val3=9/6;
        c= 'k-+';
end

for k = 1:length(fileList)
    baseFileName = fileList(k).name;
    fullFileName = fullfile(folder, baseFileName);
    params = read_chease(fullFileName);
    
    a = params.a; R0 = params.R0; B0 = params.B0;
    q = params.q; B = params.B; R = params.R; F = params.F;
    gradPsi2 = params.gradPsi2;
    r = params.r;
    nr = length(r);
    thetastar = (params.thetastar)';
    delta = params.delta(end);
    delta = params.delta(end);
    deltasTAE = [deltasTAE,params.delta(index_q_r(q,val1))];
    deltasEAE = [deltasEAE,params.delta(index_q_r(q,val2))];
    deltasNAE = [deltasNAE,params.delta(index_q_r(q,val3))];
    deltas=[deltas,delta];
    
    w = eigenmodes_chease(size,band,gradPsi2,B,B0,R,R0,F,q,thetastar,a,n,nr);
    
[TAE,reltativeTAE] = gap_at_q(val1,q,w(1,:),w(2,:));
    [EAE,reltativeEAE] = gap_at_q(val2,q,w(2,:),w(3,:));
    [NAE,reltativeNAE] = gap_at_q(val3,q,w(3,:),w(4,:));
    
    TAEgaps = [TAEgaps,reltativeTAE];
    EAEgaps = [EAEgaps,reltativeEAE];
    NAEgaps = [NAEgaps,reltativeNAE];
end


save=0; %Set 1 to save the figure, 0 otherwise

path='C:\Users\manon\Desktop\projet CSE II\figures\';

namefig = 'TAEs_ogyropsis';
figure
plot(deltasTAE,TAEgaps,c)
legend('$n=$'+ string(n)+', $q(r)=$'+ string(val1))
grid on
xlabel('$\delta$')
ylabel('Relative TAE Gap')
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

namefig = 'EAEs_ogyropsis';
figure
plot(deltasEAE,EAEgaps,c)
legend('$n=$'+ string(n)+', $q(r)=$'+ string(val2))
grid on
xlabel('$\delta$')
ylabel('Relative EAE Gap')
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

namefig = 'NAEs_ogyropsis';
figure
plot(deltasNAE,NAEgaps,c)
grid on
legend('$n=$'+string(n)+', $q(r)=$'+ string(val3))
xlabel('$\delta$')
ylabel('Relative NAE Gap')
if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end


%%

n = 3;

switch n
    case 1
        qTAE = [1.5,2.5,3.5];
        qEAE = [2,3,4];
        qNAE = [1.5,2.5,3.5];

    case 2
        qTAE = [5,7,9]./4;
        qEAE = [3,4,5]./2;
        qNAE = [5,7,9]./4;
    case 3
        qTAE = [7,9,11]./6;
        qEAE = [4,5,6]./3;
        qNAE = [7,9,11]./6;
end


TAEgaps = zeros(length(qTAE),length(fileList));
EAEgaps = zeros(length(qTAE),length(fileList));
NAEgaps = zeros(length(qTAE),length(fileList));
deltasTAE=zeros(length(qTAE),length(fileList));
deltasEAE=zeros(length(qTAE),length(fileList));
deltasNAE=zeros(length(qTAE),length(fileList));
 
for i = 1:length(qTAE)
    deltas =[];
    for k = 1:length(fileList)
    baseFileName = fileList(k).name;
    fullFileName = fullfile(folder, baseFileName);
    params = read_chease(fullFileName);
    
    a = params.a; R0 = params.R0; B0 = params.B0;
    q = params.q; B = params.B; R = params.R; F = params.F;
    gradPsi2 = params.gradPsi2;
    r = params.r;
    nr = length(r);
    thetastar = (params.thetastar)';
    delta = params.delta(end);
    deltasTAE(i,k) = params.delta(index_q_r(q,qTAE(i)));
    deltasEAE(i,k) = params.delta(index_q_r(q,qEAE(i)));
    deltasNAE(i,k) = params.delta(index_q_r(q,qNAE(i)));
    deltas=[deltas,delta];
    
    w = eigenmodes_chease(size,band,gradPsi2,B,B0,R,R0,F,q,thetastar,a,n,nr);
    
    [TAE,reltativeTAE] = gap_at_q(qTAE(i),q,w(1,:),w(2,:));
    [EAE,reltativeEAE] = gap_at_q(qEAE(i),q,w(2,:),w(3,:));
    [NAE,reltativeNAE] = gap_at_q(qNAE(i),q,w(3,:),w(4,:));
    
%     TAEgaps(i,k) = TAE;
%     EAEgaps(i,k) = EAE;
%     NAEgaps(i,k) = NAE;


    TAEgaps(i,k) = reltativeTAE;
    EAEgaps(i,k) = reltativeEAE;
    NAEgaps(i,k) = reltativeNAE;
    end
end


%%
%ERREURR
figure
plot(deltasTAE(1,:), TAEgaps(1,:),'+-')
hold on
plot(deltasTAE(2,:), TAEgaps(2,:),'+-')
hold on
plot(deltasTAE(3,:), TAEgaps(3,:),'+-')
grid on
legend('$q='+string(qTAE(1))+'$', '$q='+string(qTAE(2))+'$','$q='+string(qTAE(3))+'$')
title('$n=$'+string(n))
xlabel('$\delta$')
ylabel('Relative TAE Gap')

figure
plot(deltasEAE(1,:), EAEgaps(1,:),'+-')
hold on
plot(deltasEAE(2,:), EAEgaps(2,:),'+-')
hold on
plot(deltasEAE(3,:), EAEgaps(3,:),'+-')
grid on
legend('$q='+string(qEAE(1))+'$', '$q='+string(qEAE(2))+'$','$q='+string(qEAE(3))+'$')
title('$n=$'+string(n))
xlabel('$\delta$')
ylabel('Relative EAE Gap')

figure
plot(deltasNAE(1,:), NAEgaps(1,:),'+-')
hold on
plot(deltasNAE(2,:), NAEgaps(2,:),'+-')
hold on
plot(deltasNAE(3,:), NAEgaps(3,:),'+-')
grid on
legend('$q='+string(qNAE(1))+'$', '$q='+string(qNAE(2))+'$','$q='+string(qNAE(3))+'$')
title('$n=$'+string(n))
xlabel('$\delta$')
ylabel('Relative NAE Gap')