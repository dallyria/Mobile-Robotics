clear all
clc
syms m mo h Lo Isys d k1 k2 theta_1 theta_2 t rw b I gamma_dot theta_dot1 theta_dot2 tau_l1 tau_l2 tau_r1 tau_r2 x_dot y_dot x_dot_feedback y_dot_feedback theta_feedback theta_dot_feedback x_feedback y_feedback x_dotdot y_dotdot theta_dotdot1 theta_dotdot2;


%% Express the dynamic equations:

m=.55;
mo=0.05;
h=0.1258;
Lo=0.0500
b=0.085;
d=0.0092; %distance btw. geometric center and center of mass
%rw=0.0276;
Isys=0.0026;
k1=1/10;
save Rob_param m mo h Lo b d Isys

gamma=k1*t%+k2*t^2;
gamma_dot=diff(gamma,t);


q_dotdot=[x_dotdot;y_dotdot;theta_dotdot1;theta_dotdot2];

M=[2*(m+mo) 0  h*(m+mo)*sin(theta_1)  h*(m+mo)*sin(theta_2);...
   0 2*(m+mo) -h*(m+mo)*cos(theta_1) -h*(m+mo)*cos(theta_2);...
   h*(m+mo)*sin(theta_1) -h*(m+mo)*cos(theta_1) Isys+(m+mo)*h^2 (m+mo)*cos(theta_1-theta_2)*h^2;...
   h*(m+mo)*sin(theta_2) -h*(m+mo)*cos(theta_2) (m+mo)*cos(theta_1-theta_2)*h^2 Isys+(m+mo)*h^2];

V=[h*(m+mo)*(cos(theta_1)*theta_dot1^2+cos(theta_2)*theta_dot2^2);...
   h*(m+mo)*(sin(theta_1)*theta_dot1^2+sin(theta_2)*theta_dot2^2);...
   (m+mo)*sin(theta_1-theta_2)*(h*theta_dot2)^2;...
   -(m+mo)*sin(theta_1-theta_2)*(h*theta_dot1)^2];

Bq=[cos(theta_1) cos(theta_1)   cos(theta_2) cos(theta_2);...
    sin(theta_1) sin(theta_1)   sin(theta_2) sin(theta_2);...
    -b b -b+h*sin(theta_1-theta_2)  b+h*sin(theta_1-theta_2);...
    0 0 -b b];
%Aq = [-sin(theta) ;cos(theta);  -d]; 
%bq = -m*(x_dot*cos(theta)+y_dot*sin(theta))*theta_dot;
% Bq1=Bq(:,1:2);
% Bq2=Bq(:,3:4);
% Bq_inv=[(Bq1.'*Bq1)*Bq1.' ;(Bq2.'*Bq2)*Bq2.']

tau=[tau_l1;tau_r1;tau_l2;tau_r2];

% M*[x_dotdot;y_dotdot;theta_dotdot]+V-B*tau
% Landa=(-inv(Aq*inv(M)*Aq.')*(bq-Aq*inv(M)*Bq*tau));
qdotdot=inv(M)*(Bq*tau-V);
 
% tau_dyn=inv(Bq.'*Bq)*Bq.'*(M*[x_dotdot;y_dotdot;theta_dotdot]+V);
qdotdotfun=matlabFunction(qdotdot,'File','Qdotdot'); %Save tau to a .m file in the root folder

tau_computed=inv(Bq)*(M*[x_dotdot;y_dotdot;theta_dotdot1;theta_dotdot2]+V);
% tau_compued=Bq_inv*(M*[x_dotdot;y_dotdot;theta_dotdot1;theta_dotdot2]+V);

tau_computed_function=matlabFunction(tau_computed,'File','Tau_computed'); %Save tau to a .m file in the root folder
B_inv=matlabFunction(inv(Bq),'File','B_inverse');



save qdotdot_file qdotdot
%%
load RobotPaths2
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
q_desired=[r_t;theta_t;theta_t+pi/3];

matlabFunction(q_desired,'File','q_desired'); %Save tau to a .m file in the root folder
% matlabFunction(dr_t,'File','dR_t'); %Save tau to a .m file in the root folder
