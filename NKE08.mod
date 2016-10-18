%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo NKE basico loglinealizado. 

var c inom pic lab w y mc a g z;
varexo e_a e_g e_z;

parameters beta sigma psi phi phi_pic phi_y rho_a rho_g rho_z;

beta    = 0.995;
sigma   = 1;
psi     = 1;
phi     = 0.75;
phi_pic = 1.5;
phi_y   = 0.125;

rho_a   = 0.90;
rho_g   = 0.80;
rho_z   = 0.70;


model (linear);
c   = -1/sigma*(inom-pic(+1)) + c(+1) + 1/sigma*(1-rho_g)*g;
w   = psi*lab + (1/sigma)*c;
y   = a + lab;
pic = beta*pic(+1) + (1-phi)*(1-phi*beta)/phi*mc;
mc  = w - a;
inom= phi_pic*pic + phi_y*y + z;
y   = c;
a   = rho_a*a(-1) + e_a;
g   = rho_g*g(-1) + e_g;
z   = rho_z*z(-1) + e_z;
end;

initval;
c       = 0;
inom    = 0;
pic     = 0;
g       = 0;
w       = 0;
lab     = 0;
y       = 0;
a       = 0;
mc      = 0;
z       = 0;
end;

resid;

shocks;
var e_a; stderr 0.70;
var e_g; stderr 0.50;
var e_z; stderr 0.25;
end;

steady;

check;

stoch_simul(order=1,irf=20);
