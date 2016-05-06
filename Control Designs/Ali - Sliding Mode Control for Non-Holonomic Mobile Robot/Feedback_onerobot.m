clear all
clc
syms m d k1 k2 theta t rw b I gamma_dot theta_dot tau_l tau_r x_dot y_dot x_dot_feedback y_dot_feedback theta_feedback theta_dot_feedback x_feedback y_feedback x_dotdot y_dotdot theta_dotdot




%% Express the dynamic equations:

m=.5;
b=0.085;
d=0.03;
rw=0.0276;
I=.01;
k1=1/10;
save Rob_param m b d rw I k1

gamma=k1*t%+k2*t^2;
gamma_dot=diff(gamma,t);


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

tau=[tau_l;tau_r];

% M*[x_dotdot;y_dotdot;theta_dotdot]+V-B*tau
% Landa=(-inv(Aq*inv(M)*Aq.')*(bq-Aq*inv(M)*Bq*tau));
qdotdot=inv(M)*(Bq*tau-V-Aq*bq);
 %qdotdot=inv(M)*(Bq*[tau_l;tau_r]-V)
 
% tau_dyn=inv(Bq.'*Bq)*Bq.'*(M*[x_dotdot;y_dotdot;theta_dotdot]+V);
qdotdotfun=matlabFunction(qdotdot,'File','Qdotdot'); %Save tau to a .m file in the root folder

tau_compued=inv(Bq.'*Bq)*Bq.'*(M*[x_dotdot;y_dotdot;theta_dotdot]+V+Aq*bq);

tau_computer_function=matlabFunction(tau_compued,'File','Tau_compued'); %Save tau to a .m file in the root folder
B_inv=matlabFunction(inv(Bq.'*Bq)*Bq.','File','B_inverse');



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

matlabFunction(q_desired,'File','q_desired'); %Save tau to a .m file in the root folder
% matlabFunction(dr_t,'File','dR_t'); %Save tau to a .m file in the root folder
