%----------------------------------------------------------------------------------- 
% This code generates spline paths for mobile robot enforcing kinematic constraints
%
% This code is developed by Mohammadali Shahriari, mshahria@uoguelph.ca
%
% PhD Adviser: Dr. Mohammad Biglarbegian
% AMRCLab.com
% University of Guelph, 
% Guelph, ON, Canada
% October 2015
% -----------------------------------------------------------------------------------

%% Definition of SPline:
global Points gamma_dot Omega_max r

gamma_dot=1/10; %1/tf where tf is the travel time (sec) of the robot to reach its target.
Omega_max = 24.11; % maximum rotational speed of the wheels (rad/sec) omega*r= V,   V= 4.00 meters / 6 seconds 
r= 0.0276 ; % Wheel's radius (m)

syms T R g  g_dot real
ndegree=3; % order of splines used for fitting to robot's path.

A = sym('a',[2,ndegree+1]); % polynomial variables of x and y
assume(A,'real'); % assumptions;
gv=g.^(0:ndegree); % SPlines are in terms of \gamma in [0,1] ,  \gamma= k*t   0<k<= 1/t_final

ri=A*[gv].'; % position of the robot 
v=g_dot*diff(ri,g); % v= [x_dot; y_dot]
Vn=sqrt(v.'*v); 
dV=diff(Vn,g);
theta=atan2(v(2),v(1));
omega=g_dot*diff(theta,g);

% Save Path as a function of SPline coefficients to matlab files
Ri=matlabFunction(ri,'vars',{[A(1,:),A(2,:)],g},'file','Ri');
Vi=matlabFunction(v,'vars',{[A(1,:),A(2,:)],g,g_dot},'file','Vi');
dV=matlabFunction(dV,'vars',{[A(1,:),A(2,:)],g,g_dot},'file','dV');
absVi=matlabFunction(Vn,'vars',{[A(1,:),A(2,:)],g,g_dot},'file','absVi');
Omg=matlabFunction(omega,'vars',{[A(1,:),A(2,:)],g,g_dot},'file','Omegafun');
%% Finding SPline coefficients with Wheeled mobile robot constraints


%====Rout 1=====================
for i=1:5
Poi{i}=[.20 1.15 0;...
    1.70 .75-i*0.1 0.3;...
     2.70 1.0 .7;...
    3.80 1.00 1];
end


for iPo=1:size(Poi,2)
Points=[];    
Points=Poi{iPo};


options=optimset('Display','iter-detailed')%'UseParallel',true);
[Aspl, fval, exitflag]=fmincon(@SPline_Objective,(zeros(1,2*ndegree+2)),[],[],[],[],[],[],@WMR_constraints,options)

p_x=Aspl(1,1:ndegree+1);
p_y=Aspl(1,ndegree+2:end);
gamma_array=linspace(0,1,100);
Rbt_pos=Ri([p_x,p_y],gamma_array);
Pol{iPo}=[p_x;p_y];
end

axis tight
axis equal
title Eenvironment
xlabel x(m)
ylabel y(m)
savefig('Environment.fig')


