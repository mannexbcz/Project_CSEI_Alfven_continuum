% -------------------------------------------------------------------------
% This script compares the values of R and Z obtained with and without
% an equidistant mesh in Thetastar
% -------------------------------------------------------------------------

%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;
%%
filename = '';

psichease   = hdf5read(filename,'/data/grid/PSI');
psitilde    = psichease/psichease(length(psichease));
schease     = sqrt(psichease/psichease(length(psichease)));
chichease   = hdf5read(filename,'/data/grid/CHI');

%% R,Z from CHEASE

R = hdf5read(filename,'/data/var2d/R');
Z = hdf5read(filename,'/data/var2d/Z');
Ra = R(end,:);
Za = Z(end,:);

figure
plot(Ra,Za,'k')
axis equal
xlabel('$R$')
ylabel('$Z$')
%%
Ra0 = (1/2)*(min(Ra)+max(Ra));
%% R,Z Geometric
params = read_chease(filename);

Rgeom   = params.Rgeom;
R0      = params.R0;
Delta   = params.Delta;
r       = params.r;
a       = params.a;
delta   = params.delta;
kappa   = params.kappa;
q       = params.q;
B0      = params.B0;

theta = linspace(0,2*pi,512);

%% déterminer locations of q(r)={1.5,2,2.5}
q = hdf5read(filename,'/data/var1d/q');
figure 
plot(q,'k')
hold on
plot(1.5.*ones(1,length(r)),'b-')
hold on
plot(2.*ones(1,length(r)),'b-')
hold on
plot(2.5.*ones(1,length(r)),'b-')
xlabel('r')
ylabel('q(r)')

% q(r)=1.5 => r=0.154646, index 175 // redone: 175
% q(r)=2   => r=0.191345, index 206 // redone: 206
% q(2)=2.5 => r=0.214323, index 223 // redone: 222
%%
% Plasma Boundary
Rg = R0 + Delta(end)+ a.*cos(theta+asin(delta(end)).*sin(theta));
Zg = kappa(end).*a.*sin(theta);
% q = 1.5
R15g =R0 + Delta(175)+ r(175).*cos(theta+asin(delta(175)).*sin(theta));
Z15g = kappa(175).*r(175).*sin(theta);
R15a = R(175,:);
Z15a = Z(175,:);
% q = 2
R2g =R0 + Delta(206)+ r(206).*cos(theta+asin(delta(206)).*sin(theta));
Z2g = kappa(206).*r(206).*sin(theta);
R2a = R(206,:);
Z2a = Z(206,:);
% q = 2.5
R25g =R0 + Delta(223)+ r(223).*cos(theta+asin(delta(223)).*sin(theta));
Z25g = kappa(223).*r(223).*sin(theta);
R25a = R(223,:);
Z25a = Z(223,:);
%%
save = 0;
namefig = 'RZsurfaces';
path='C:\Users\manon\Desktop\projet CSE II\figures\';

figure
plot(Ra,Za,'-k')
hold on 
plot(Rg,Zg,'b-')
hold on
plot(R15a,Z15a,'-k')
hold on 
plot(R15g,Z15g,'b-')
hold on
plot(R2a,Z2a,'-k')
hold on 
plot(R2g,Z2g,'b-')
hold on
plot(R25a,Z25a,'-k')
hold on 
plot(R25g,Z25g,'b-')
axis equal
xlabel('$R$')
ylabel('$Z$')
legend('Chease','Geometric')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
 
%% R,Z wrt thetastar
npoints = 512;
j=0;
thetastar = zeros(4,length(theta));
for i=[175,206,223,length(r)]
    j=j+1;
    epsilon=params.r(i)/params.R0;
    delta = params.Delta(i);
    k = params.kappa(i);
    d = params.delta(i);
    kprim = params.kprime(i);
    deltaprim = params.deltaprime(i);
    dprim = params.dprime(i);
    qbar = qbar_triang(params.q(i),params.r(i),k,kprim,delta,deltaprim,d,dprim,params.R0,npoints);
    
    thetastar(j,:) = theta_star_triang(r(i),theta,kprim,k,delta,deltaprim,d,dprim,q(i),qbar,B0,R0,npoints); 
end

%% Plot R wrt thetastar

figure
C15 = plot(chichease,R15a, 'k-');
hold on
G15 = plot(thetastar(1,:),R15g, 'k--');
hold on
C2 = plot(chichease,R2a, 'k-');
hold on
G2 = plot(thetastar(2,:),R2g, 'k--');
hold on
C25 = plot(chichease,R25a, 'k-');
hold on
G25 = plot(thetastar(3,:),R25g, 'k--');
xlabel('$\theta^*$')
ylabel('$R$')
legend('CHEASE','Geometric')

%% Plot Z wrt thetastar
figure
C15 = plot(chichease,Z15a, 'k-');
hold on
G15 = plot(thetastar(1,:),Z15g, 'k--');
hold on
C2 = plot(chichease,Z2a, 'k-');
hold on
G2 = plot(thetastar(2,:),Z2g, 'k--');
hold on
C25 = plot(chichease,Z25a, 'k-');
hold on
G25 = plot(thetastar(3,:),Z25g, 'k--');
xlabel('$\theta^*$')
ylabel('$Z$')
legend('CHEASE','Geometric')