clear all
clc
syms m d k1 k2 theta t rw b I gamma_dot theta_dot tau_l tau_r x_dot y_dot x_dot_feedback y_dot_feedback theta_feedback theta_dot_feedback x_feedback y_feedback x_dotdot y_dotdot theta_dotdot
syms x_desired theta_desired x_dot_desired y_desired y_dot_desired theta_dot_desired  x y lambda_1 lambda_2 x_dotdot_desired y_dotdot_desired theta_dotdot_desired
 



%% Express the dynamic equations:

m=.5;
b=0.085;
d=0.03;
rw=0.0276;
I=.01;
k1=1/10;
save Rob_param m b d rw I k1

gamma=k1*t;%+k2*t^2;
gamma_dot=diff(gamma,t);


vc=sqrt(x_dot^2+y_dot^2);
wc=theta_dot;



vc_dot=sqrt(x_dotdot^2+y_dotdot^2);
wc_dot=theta_dotdot;

vc_dot_desired=sqrt(x_dotdot_desired^2+y_dotdot_desired^2)
wc_dot_desired=theta_dotdot_desired


z=[vc;wc];
zdot=[vc_dot;wc_dot];

M=[m 0  d*m*sin(theta);...
    0 m -d*m*cos(theta);...
    d*m*sin(theta) -d*m*cos(theta) I];


V=[d*m*cos(theta)*theta_dot^2;...
   d*m*sin(theta)*theta_dot^2;...
   0];

Bq=1/rw *[cos(theta) cos(theta);...
    sin(theta) sin(theta);...
    b -b];
Aq = [-sin(theta) ;cos(theta);  -d]; 
bq = -m*(x_dot*cos(theta)+y_dot*sin(theta))*theta_dot;

S=[cos(theta) 0; sin(theta) 0 ; 0  1];

Hq= inv(S.'*Bq)*S.'*M*S;
S_dot=[-theta_dot*sin(theta) 0;theta_dot* cos(theta) 0 ; 0  0];

Gq=inv(S.'*Bq)*S.'*(M*S_dot*z+V)

tau=[tau_l;tau_r];

% M*[x_dotdot;y_dotdot;theta_dotdot]+V-B*tau
% Landa=(-inv(Aq*inv(M)*Aq.')*(bq-Aq*inv(M)*Bq*tau));
qdotdot=inv(M)*(Bq*tau-V-Aq*bq);
 %qdotdot=inv(M)*(Bq*[tau_l;tau_r]-V)
 
% tau_dyn=inv(Bq.'*Bq)*Bq.'*(M*[x_dotdot;y_dotdot;theta_dotdot]+V);
qdotdotfun=matlabFunction(qdotdot,'File','Qdotdot'); %Save tau to a .m file in the root folder

tau_fdbckln=Hq*[vc_dot_desired;wc_dot_desired]+Gq;
matlabFunction(tau_fdbckln,'File','Tau_fdbckln')
matlabFunction(Hq,'File','Hqfun');


%Sliding 
rho_desired=sqrt(x_desired^2+y_desired^2);
rhoc=sqrt(x^2+y^2)
rho_e=rho_desired-rhoc;
rho_dot_e=sqrt(x_dot_desired^2+y_dot_desired^2)-sqrt(x_dot^2+y_dot^2);
phi_c=atan2(y,x);
phi_desired=atan2(y_desired,x_desired);
s1=rho_dot_e+lambda_1 *rho_e;
s2=(theta_dot_desired-theta_dot)+lambda_2*(theta_desired-theta)+sign(theta_desired-theta)*abs(phi_desired-phi_c);

syms Q1 P1 Q2 P2

vc_desired=sqrt(x_dot_desired^2+y_dot_desired^2);
wc_desired=theta_dot_desired;

Rqz=vc^2*sin(phi_c-theta)^2/rhoc+vc*wc*sin(phi_c-theta);
Rqz_desired=vc_desired^2*sin(phi_desired-theta_desired)^2/rho_desired+vc_desired*wc_desired*sin(phi_desired-theta_desired);


u1=(1/cos(phi_c-theta))*(-Q1*s1-P1*sign(s1)-lambda_1*rho_dot_e+vc_dot_desired*cos(phi_desired-theta_desired)+Rqz_desired-Rqz)-vc_dot_desired;
u2=-Q2*s2-P2*sign(s2)-lambda_2 *(theta_dot_desired-theta_dot)-sign(theta_desired-theta)*abs(atan2(y_dot_desired,x_dot_desired)-atan2(y_dot,x_dot));

u=[u1;u2];

matlabFunction(u,'File','U_sliding')

save qdotdot_file qdotdot
%%
load RobotPaths
Poly=Pol{1,1};
Pol_A=Poly(1,:)
Pol_B=Poly(2,:)

D=length(Pol_A)-1;%order of spline
A=sym('a',[1,D+1]);
B=sym('b',[1,D+1]);
r=[A;B]*[gamma.^(0:D)].';
r_t=[Pol_A;Pol_B]*[gamma.^(0:D)].';
dr_t=diff(r_t,t);
theta_t=atan2(dr_t(2),dr_t(1));
q_desired=[r_t;theta_t];
q_dotdot_desired=diff(diff(q_desired,t),t);
matlabFunction(diff(q_desired,t),'File','q_dot_desired'); 

matlabFunction(q_desired,'File','q_desired'); 
matlabFunction(q_dotdot_desired,'File','q_dotdot_desired'); 
%Save tau to a .m file in the root folder
% matlabFunction(dr_t,'File','dR_t'); %Save tau to a .m file in the root folder
