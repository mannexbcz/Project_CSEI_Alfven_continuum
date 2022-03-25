% filename=input('Enter CHEASE h5 output file name: ','s')
filename='ogyropsi';
psichease=hdf5read(filename,'/data/grid/PSI');
psitilde = psichease/psichease(length(psichease));
schease=sqrt(psichease/psichease(length(psichease)));
chichease=hdf5read(filename,'/data/grid/CHI');
npsi=length(psichease);
nchi=length(chichease);

q  = hdf5read(filename,'/data/var1d/q');
F  = hdf5read(filename,'/data/var1d/f');
rho_tor  = hdf5read(filename,'/data/var1d/rho_tor');
ageom    = hdf5read(filename,'/data/var1d/ageom');
Rgeom    = hdf5read(filename,'/data/var1d/Rgeom');
kappa    = hdf5read(filename,'/data/var1d/kappa');
delta_lower=hdf5read(filename,'/data/var1d/delta_lower');
delta_upper=hdf5read(filename,'/data/var1d/delta_upper');
amid = ageom(length(ageom));
Rmid = Rgeom(1);

R        = hdf5read(filename,'/data/var2d/R');
Z        = hdf5read(filename,'/data/var2d/Z');
B        = hdf5read(filename,'/data/var2d/B');
Jacobian = hdf5read(filename,'/data/var2d/Jacobian');
dPsidR   = hdf5read(filename,'/data/var2d/dPsidR');
dPsidZ   = hdf5read(filename,'/data/var2d/dPsidZ');

gradPsi2 = dPsidR.^2+dPsidZ.^2;

% we chose to normalize 
%   R0 = R_magnetic axis
%   B0 = B_magnetic axis
R0=Rgeom(1);
B0=F(1)/R0;

aeq=(B0^2/R0^2)*gradPsi2.*R.^2./B.^2;
beq=gradPsi2./(B.^2.*R.^2);








