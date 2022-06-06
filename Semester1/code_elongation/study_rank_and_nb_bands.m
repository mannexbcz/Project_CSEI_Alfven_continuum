% -------------------------------------------------------------------------
% This script studies the impact on the convergence of the rank and the
% number of bands of the matrices M and N considered, in the case of an
% elongated and shifted equilibriuum, and in the first order approximation.
% The computation time with respect to these two parameters is also computed.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
B0=1;       % magnetic field
F = R0*B0;  % toroidal flux current
nr = 1000;  % number of magnetic surfaces considered
npoints=100;% number of mesh points used for the integration
q = @(r) 1+2*(r/a).^2; % safety factor

% Elongation
k_a=1;
kprime=0;
kfct = @(r) kprime.*(r-a) + k_a;

%Shift
deltaprime= 0;
deltafct = @(r)  -deltaprime*a + deltaprime *r;

%% Convergence with respect to the number of bands of the matrices
size=10;     % rank of the matrices M and N
bands=[3,5,7,9,11,13,15,17,19];  % number of bands of the matrices M and N

gaps=zeros(1,length(bands));
time_bands = zeros(1,length(bands));

r=linspace(0,a,nr);

for s=1:length(bands)

band=bands(s);
    
w=zeros(size,length(r));
w2=zeros(size,length(r));
time = 0;
for i=1:length(r) 
    tic

    epsilon=r(i)/R0;
    delta = deltafct(r(i));
    k = kfct(r(i));
    qbar = qbar_ordre1(q(r(i)),r(i),k,kprime);

    [M,N]=matrices_ordre1(r(i),epsilon,k,kprime,delta,deltaprime,q(r(i)),qbar,R0,B0,n,size,band,npoints);
    
    w2(:,i)=eig(N,M); 
    for j=1:size 
        if w2(j,i)<0
            w2(j,i)=nan;
        end
    end
    w2(:,i)=sort(w2(:,i));
    w(:,i)=sqrt(w2(:,i));
    
    time = time + toc;
end

[pks1,locs1] = findpeaks(real(w(1,:)),'MinPeakDistance',10);
[pks2,locs2] = findpeaks(-real(w(2,:)),'MinPeakDistance',10);

gaps(s)=w(2,locs2(1))-w(1,locs1(1));
time_bands(s)=time/length(r);
end

%% Figure

save=1; %Set 1 to save the figure, 0 otherwise
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
%% Figure Time 

save=1; %Set 1 to save the figure, 0 otherwise
namefig=['time_bands_ordre1'];
path='C:\Users\manon\Desktop\projet CSE I\figures\elongation\';

figure
plot(bands,time_bands,'k +')
grid on
xlabel('Bands')
ylabel({'Computation time';'per magnetic surface [s]'})

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end

%% Convergence with respect to the rank of the matrices

sizes=[2:40];  % rank of the matrices M and N
band = 5; % number of bands of the matrices M and N

gaps=zeros(1,length(sizes));
time_size = zeros(1,length(sizes));

for s=1:length(sizes)

size=sizes(s);
    
w=zeros(size,length(r));
w2=zeros(size,length(r));
time = 0;
for i=1:length(r) 
    tic

    epsilon=r(i)/R0;
    delta = deltafct(r(i));
    k = kfct(r(i));
    qbar = qbar_ordre1(q(r(i)),r(i),k,kprime);

    [M,N]=matrices_ordre1(r(i),epsilon,k,kprime,delta,deltaprime,q(r(i)),qbar,R0,B0,n,size,band,npoints);

    w2(:,i)=eig(N,M); 
    for j=1:size 
        if w2(j,i)<0
            w2(j,i)=nan;
        end
    end
    w2(:,i)=sort(w2(:,i));
    w(:,i)=sqrt(w2(:,i));
    time = time +toc;
end

[pks1,locs1] = findpeaks(real(w(1,:)),'MinPeakDistance',10);
[pks2,locs2] = findpeaks(-real(w(2,:)),'MinPeakDistance',10);

gaps(s)=w(2,locs2(1))-w(1,locs1(1));
time_size(s) = time/length(r) ;
end

%% Reference

size=100;
band = 5;

w=zeros(size,length(r));
w2=zeros(size,length(r));
for i=1:length(r)  
    epsilon=r(i)/R0;
    delta = deltafct(r(i));
    k = kfct(r(i));
    qbar = qbar_ordre1(q(r(i)),r(i),k,kprime);

    [M,N]=matrices_ordre1(r(i),epsilon,k,kprime,delta,deltaprime,q(r(i)),qbar,R0,B0,n,size,band,npoints);

    w2(:,i)=eig(N,M); 
    for j=1:size 
        if w2(j,i)<0
            w2(j,i)=nan;
        end
    end
    w2(:,i)=sort(w2(:,i));
    w(:,i)=sqrt(w2(:,i));
end

[pks1,locs1] = findpeaks(real(w(1,:)),'MinPeakDistance',10);
[pks2,locs2] = findpeaks(-real(w(2,:)),'MinPeakDistance',10);

gaps_ref=w(2,locs2(1))-w(1,locs1(1));

%% Figure

save=1; %Set 1 to save the figure, 0 otherwise
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

%% Figure Time 

save=1; %Set 1 to save the figure, 0 otherwise
namefig=['time_size_ordre1'];
path='C:\Users\manon\Desktop\projet CSE I\figures\elongation\';

figure
plot(sizes,time_size,'k +')
grid on
xlabel('Rank')
ylabel({'Computation time';'per magnetic surface [s]'})

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end