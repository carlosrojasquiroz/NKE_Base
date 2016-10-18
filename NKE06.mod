%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo NKE basico loglinealizado. 
// El Estado Estacionario es cero debido a la loglinealización manual.
// Modelo semiestructural con componente backward looking
// (c) Carlos Rojas Quiroz 

var ygap inom pic g u z;
varexo ghat uhat zhat;
parameters phi sigma psi beta lambda Psi gamma_x alpha_pic rho_i phi_x phi_pic rho mu theta;

phi         =0.60;
sigma       =1.00;
psi         =0.75; 
beta        =0.99;
rho         =0.75;
mu          =0.75;
theta       =0.75;
lambda      =phi*psi*(1-(1-psi)*beta)/(1-psi);
Psi         =1/sigma;
gamma_x     =0.45;
alpha_pic   =0.45;
rho_i       =0.70;
phi_x       =0.25;
phi_pic     =1.50;

model(linear);
ygap    =gamma_x*ygap(-1) + (1-gamma_x)*ygap(+1)-Psi*(inom-pic(+1))+g;
pic     =alpha_pic*pic(-1)+ (1-alpha_pic)*pic(+1)+lambda*ygap+u;
inom    =rho_i*inom(-1)+(1-rho_i)*(phi_x*ygap+phi_pic*pic)+z;
g       =mu*g(-1)+ghat;
u       =rho*u(-1)+uhat;
z       =theta*z(-1)+zhat;
end;

shocks;
var ghat; stderr 0.01;
var uhat; stderr 0.01;
var zhat; stderr 0.01;
end;


stoch_simul(irf=20, nograph);
