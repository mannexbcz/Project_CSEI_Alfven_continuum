%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;
%%
filename = 'ogyropsi_positive_triangularity.h5';

psichease   = hdf5read(filename,'/data/grid/PSI');
psitilde    = psichease/psichease(length(psichease));
schease     = sqrt(psichease/psichease(length(psichease)));
chichease   = hdf5read(filename,'/data/grid/CHI');

%% R,Z from CHEASE

R = hdf5read(filename,'/data/var2d/R');
Z = hdf5read(filename,'/data/var2d/Z');
Ra = R(end,:);
Za = Z(end,:);

Ra0 = (1/2)*(min(Ra)+max(Ra));
%% R,Z Geometric

Rgeom   = hdf5read(filename,'/data/var1d/Rgeom');
R0      =  Rgeom(end);
Delta   = Rgeom - Rgeom(end);
Volume  = hdf5read(filename,'/data/var1d/Volume');
rho_vol = sqrt(Volume/Volume(length(Volume)));
r = max(R,[],2)-Rgeom;
a       = r(end);
delta_lower      = hdf5read(filename,'/data/var1d/delta_lower');
delta_upper      = hdf5read(filename,'/data/var1d/delta_upper');
delta            = (1/2)*(delta_lower+delta_upper);
kappa       = hdf5read(filename,'/data/var1d/kappa');

theta = linspace(0,2*pi,100);

%% déterminer les r tq q(r)=1.5 etc pour ploter R,Z au nv des gaps
q     = hdf5read(filename,'/data/var1d/q');
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

% q(r)=1.5 => r=0.154646, index 175
% q(r)=2   => r=0.191345, index 206
% q(2)=2.5 => r=0.214323, index 223
%%
% Dernière surface
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
%%
% params = read_chease(filename);
% 
% Rfctr = R0 + Delta+ r.*cos(chichease(36)+asin(delta).*sin(chichease(36)));
% 
% figure
% plot(r,R(:,36))
% hold on
% plot(r,Rfctr)