%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo NKE basico loglinealizado. 

var c inom pic lab w y mc m pm z;
varexo e_pm e_z;

parameters beta sigma eta phi rho_i phi_pic phi_y rho_pm rho_z;

beta    = 0.99;
sigma   = 1;
eta     = 1;
phi     = 0.75;
rho_i   = 0.70;
phi_pic = 1.50;
phi_y   = 0.25;
rho_z   = 0.75;
rho_pm  = 0.50;

model (linear);
c   = -1/sigma*(inom-pic(+1)) + c(+1);
w   = eta*lab + (1/sigma)*c;
y   = z + lab;
pic = beta*pic(+1) + (1-phi)*(1-phi*beta)/phi*mc;
mc  = w - z;
inom= rho_i*inom(-1)+(1-rho_i)*(phi_pic*pic + phi_y*y) + pm;
m   =sigma*c-beta*inom;
y   = c;
z   = rho_z*z(-1) + e_z;
pm  = rho_pm*pm(-1) + e_pm;
end;

initval;
c       = 0;
inom    = 0;
pic     = 0;
w       = 0;
lab     = 0;
y       = 0;
pm      = 0;
mc      = 0;
m       = 0;
z       = 0;
end;

resid;

shocks;
var e_pm; stderr 0.01;
var e_z; stderr 0.01;
end;

steady;

check;

stoch_simul(order=1, nograph);
