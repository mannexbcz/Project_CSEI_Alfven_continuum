function parameters = read_chease(filename)
% Paramètres figures
% [a,R0,B0,k_a,kprime,deltaprime,dprime,q]
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%  filename = 'ogyropsi_negative_triangularity.h5';
fig = 0;

psichease   = hdf5read(filename,'/data/grid/PSI');
psitilde    = psichease/psichease(length(psichease));
schease     = sqrt(psichease/psichease(length(psichease)));
chichease   = hdf5read(filename,'/data/grid/CHI');
parameters.thetastar = chichease;
npsi        = length(psichease);
nchi        = length(chichease);

R           = hdf5read(filename,'/data/var2d/R');
parameters.R = R;
Rgeom   = hdf5read(filename,'/data/var1d/Rgeom');
R0      =  Rgeom(end);
r = max(R,[],2)-Rgeom;
parameters.r = r;
parameters.a = r(end);

%% Computation of rho_vol
% Volume      = hdf5read(filename,'/data/var1d/Volume');
% rho_vol = sqrt(Volume/Volume(length(Volume)));
% parameters.rho_vol     = rho_vol;

%% Elongation
kappa = hdf5read(filename,'/data/var1d/kappa');
parameters.kappa = kappa;
% fit = polyfit(r,kappa,3);
% parameters.kprime     = polyval(polyder(fit),r);
% parameters.k_a         = polyval(fit, r(end));

parameters.kprime = zeros(length(r),1);
parameters.kprime(2:end-1,1)= (kappa(3:end)-kappa(1:end-2))./(r(3:end)-r(1:end-2));
parameters.kprime(1,1)= (-3*kappa(1)+4*kappa(2)-kappa(3))/(r(3)-r(1));
parameters.kprime(end,1)=(kappa(end-2)-4*kappa(end-1)+3*kappa(end))/(r(end)-r(end-2));
% parameters.k_a = kappa(end);

if fig
    figure
    plot(r, kappa, 'k+')
    hold on
    plot(r, parameters.kprime.*(r-parameters.a)+parameters.k_a)
%     hold on
%     plot(r, kprimetest.*(r-parameters.a)+parameters.k_a,'r')
    xlabel('$r$')
    ylabel('$\kappa$')
end
%% Shift
% Rmin    = min(R,[],2);
% Rmax    = max(R,[],2);
% Rgeom   = (1/2)*(Rmin+Rmax);
Delta = Rgeom - Rgeom(end);
parameters.Delta = Delta;
% fit = polyfit(r,Delta,3);
% parameters.deltaprime = polyval(polyder(fit),r);

parameters.deltaprime = zeros(length(r),1);
parameters.deltaprime(2:end-1,1)= (Delta(3:end)-Delta(1:end-2))./(r(3:end)-r(1:end-2));
parameters.deltaprime(1,1)= (-3*Delta(1)+4*Delta(2)-Delta(3))/(r(3)-r(1));
parameters.deltaprime(end,1)=(Delta(end-2)-4*Delta(end-1)+3*Delta(end))/(r(end)-r(end-2));

if fig
    figure
    plot(r, Delta, 'k+')
    hold on
    plot(r, parameters.deltaprime.*(r-parameters.a))
    xlabel('$r$')
    ylabel('$\Delta$')
end

%% Triangularity
delta_lower      = hdf5read(filename,'/data/var1d/delta_lower');
delta_upper      = hdf5read(filename,'/data/var1d/delta_upper');
delta = (1/2)*(delta_lower+delta_upper);
parameters.delta = delta;
% fit = polyfit(r,delta,3);
% parameters.dprime = polyval(polyder(fit),r);
parameters.dprime = zeros(length(r),1);
parameters.dprime(2:end-1,1)= (delta(3:end)-delta(1:end-2))./(r(3:end)-r(1:end-2));
parameters.dprime(1,1)= (-3*delta(1)+4*delta(2)-delta(3))/(r(3)-r(1));
parameters.dprime(end,1)=(delta(end-2)-4*delta(end-1)+3*delta(end))/(r(end)-r(end-2));

if fig
    figure
    plot(r, delta, 'k+')
    hold on
    plot(r, parameters.dprime.*r)
    xlabel('$r$')
    ylabel('$\delta$')
end
%%
parameters.q     = hdf5read(filename,'/data/var1d/q');
parameters.F     = hdf5read(filename,'/data/var1d/f');
parameters.R0    = Rgeom(1);
parameters.B0    = parameters.F(1)/parameters.R0;
parameters.B     = hdf5read(filename,'/data/var2d/B');

dPsidR = hdf5read(filename,'/data/var2d/dPsidR');
dPsidZ = hdf5read(filename,'/data/var2d/dPsidZ');

parameters.gradPsi2 = dPsidR.^2+dPsidZ.^2;
 return