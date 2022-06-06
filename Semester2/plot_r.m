% For a given equilibrium, this scripts plots the values of rho_vol,
% rho_tor and s_chease wrt r.
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%%

[filename, pathname] = uigetfile( {'*.*','All Files (*.*)'},'Pick the CHEASE equilibria file');
if(filename==0)
    warning('No file selected.')
end
path = append(pathname,filename);
params = read_chease(path);

%%

figure
plot(params.r./params.a, params.rho_vol)
hold on
plot(params.r./params.a, params.rho_tor./params.rho_tor(end))
hold on
plot(params.r./params.a, params.schease)
hold on
plot(params.r./params.a,params.r./params.a,'k--')
xlabel('$r/a$')
grid on
legend('$\rho_\mathrm{vol}$','$\rho_\mathrm{tor}$','$s_\mathrm{chease}$','$y=x$')