%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo NKE basico en niveles. 
// (c) Carlos Rojas Quiroz 

var z pic inom picS mc w nuP lab y c x1 x2 m;
varexo e_z;
parameters sigma beta psi eta theta phi epsilon rho_z rho_i phi_y phi_pic
z_ss pic_ss i_ss picS_ss mc_ss w_ss nuP_ss n_ss y_ss c_ss x1_ss x2_ss m_ss;
sigma   =1;
beta    =0.99;
psi     =1;
eta     =1;
theta   =1;
phi     =0.75;
epsilon =10;
rho_z   =0.75;
rho_i   =0.70;
phi_y   =0.25;
phi_pic =1.50;
z_ss    =1;
pic_ss  =0.02;
i_ss    =((1+pic_ss)/beta-1);
picS_ss =(((1+pic_ss)^(1-epsilon)-phi)/(1-phi))^(1/(1-epsilon))-1;
mc_ss   =(1+picS_ss)/(1+pic_ss)*(epsilon-1)/epsilon*(1-phi*beta*(1+pic_ss)^epsilon)/(1-phi*beta*(1+pic_ss)^(epsilon-1));
w_ss    =mc_ss*z_ss;
nuP_ss  =(1-phi)*(1+picS_ss)^(-epsilon)*(1+pic_ss)^epsilon/(1-((1+pic_ss)^epsilon)*phi);
n_ss    =((nuP_ss/z_ss)^sigma*w_ss/psi)^(1/(eta+sigma));
y_ss    =z_ss*n_ss/nuP_ss;
c_ss    =y_ss;
x1_ss   =c_ss^(1-sigma)*mc_ss/(1-phi*beta*(1+pic_ss)^epsilon);
x2_ss   =c_ss^(1-sigma)/(1-phi*beta*(1+pic_ss)^(epsilon-1));
m_ss    =theta*c_ss^sigma*(1+i_ss)/i_ss;

model;
exp(c)^(-sigma)         =beta*exp(c(+1))^(-sigma)*(1+exp(inom))/(1+exp(pic(+1)));
psi*exp(lab)^eta        =exp(c)^(-sigma)*exp(w);
exp(m)                  =theta*(1+exp(inom))/exp(inom)*exp(c)^sigma;
exp(mc)                 =exp(w)/exp(z);
exp(c)                  =exp(y);
exp(y)                  =exp(z)*exp(lab)/exp(nuP);
exp(nuP)                =(1-phi)*(1+exp(picS))^(-epsilon)*(1+exp(pic))^epsilon+(1+exp(pic))^epsilon*phi*exp(nuP(-1));
(1+exp(pic))^(1-epsilon)=(1-phi)*(1+exp(picS))^(1-epsilon)+phi;
1+exp(picS)             =epsilon/(epsilon-1)*(1+exp(pic))*exp(x1)/exp(x2);
exp(x1)                 =exp(c)^(-sigma)*exp(mc)*exp(y)+phi*beta*(1+exp(pic(+1)))^epsilon*exp(x1(+1));
exp(x2)                 =exp(c)^(-sigma)*exp(y)+phi*beta*(1+exp(pic(+1)))^(epsilon-1)*exp(x2(+1));
z                       =rho_z*z(-1)+e_z;
exp(inom)/i_ss          =(exp(inom(-1))/i_ss)^rho_i*((exp(y)/y_ss)^phi_y*(exp(pic)/pic_ss)^phi_pic)^(1-rho_i);

end;

steady_state_model;
z       =log(z_ss);
pic     =log(pic_ss);
inom    =log(i_ss);
picS    =log(picS_ss);
mc      =log(mc_ss);
w       =log(w_ss); 
nuP     =log(nuP_ss);
lab     =log(n_ss);
y       =log(y_ss);
c       =log(c_ss);
x1      =log(x1_ss);
x2      =log(x2_ss);
m       =log(m_ss);
end;

shocks;
var e_z; stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, nograph);
