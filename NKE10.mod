var c lab w y mcH pic pih rer inom z mon g ph_p pih da dah daf x yf;
varexo e_z e_mon e_g e_yf;
parameters sigma eta alpha alphaH eta_c eta_star beta phi rho_i rho_yf
phi_y phi_pic rho_z rho_mon rho_g G_Y C_Y;
sigma   =1.00;
eta     =1.00;
alpha   =0.40;
alphaH  =0.65;
eta_c   =1.00;
eta_star=1.00;
beta    =0.99;
phi     =0.75;
rho_i   =0.70;
phi_y   =0.10;
phi_pic =1.50; 
rho_z   =0.95;
rho_mon =0.50; 
rho_g   =0.50;
rho_yf  =0.75;
G_Y     =0.15;
C_Y     =1-G_Y;

model(linear);
da      = C_Y*c + G_Y*g;
dah     = da - eta_c*ph_p;
daf     = da - eta_c*rer;
x       = yf - eta_star*(ph_p-rer);
0       = (1-alpha)*ph_p+alpha*rer;
c       = c(+1) - 1/sigma*(inom - pic(+1));
eta*lab = w - sigma*c;
y       = (1-alpha)*dah+alpha*x;
y       = z + lab;
mcH     = w - z;
pih     = beta*pih(+1) + (1-phi)*(1-beta*phi)/phi*(mcH-ph_p);
pic     = pih + alpha/(1-alpha)*(rer-rer(-1));
rer     = sigma*(c-yf);
inom    = phi_y*y + phi_pic*pic + mon;
z       = rho_z*z(-1) + e_z;
mon     = rho_mon*mon(-1) + e_mon;
g       = rho_g*g(-1) + e_g;
yf      = rho_yf*yf(-1) + e_yf;
end;

shocks;
var e_z;   stderr 0.01;
var e_mon; stderr 0.01;
var e_g;   stderr 0.01;
var e_yf;  stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, nograph);
