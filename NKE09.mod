var c lab inv w y r kap mc pic inom z mon g;
predetermined_variables kap;
varexo e_z e_mon e_g;
parameters sigma eta delta alpha beta phi rho_i 
phi_y phi_pic rho_z rho_mon rho_g G_Y I_Y C_Y;
sigma  =1.00;
eta    =1.00;
delta  =0.025;
alpha  =0.65;
beta   =0.99;
phi    =0.75;
rho_i  =0.70;
phi_y  =0.10;
phi_pic=1.50; 
rho_z  =0.95;
rho_mon=0.50; 
rho_g  =0.50;
G_Y    =0.15;
I_Y    =delta*(1-alpha)*beta/(1-beta+beta*delta);
C_Y    =1-G_Y-I_Y;

model(linear);
c       = c(+1) - 1/sigma*(inom - pic(+1));
eta*lab = w - sigma*c;
kap(+1) = (1-delta)*kap + delta*inv;
r(+1)   = inom - pic(+1);
w       = mc + y - lab;
y       = C_Y*c + G_Y*g + I_Y*inv;
y       = z + alpha*lab+ (1-alpha)*kap;
mc      = alpha*w + (1-alpha)*r - z;
pic     = beta*pic(+1) + (1-phi)*(1-beta*phi)/phi*mc;
inom    = rho_i*inom(-1) + (1-rho_i)*phi_y*y + (1-rho_i)*phi_pic*pic + mon;
z       = rho_z*z(-1) + e_z;
mon     = rho_mon*mon(-1) + e_mon;
g       = rho_g*g(-1) + e_g;
end;

shocks;
var e_z;   stderr 0.01;
var e_mon; stderr 0.01;
var e_g;   stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, nograph);
