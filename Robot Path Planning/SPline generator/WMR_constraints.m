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


function  [c,ceq]= WMR_constraints(x)
ceq=[];
c=[];
global gamma_dot Omega_max r

%find the extermums of |V|
syms g
DV=dV(x,g,gamma_dot);
Sol=vpasolve(DV,g,[0 1]);

if ~isempty(Sol)
    Crit_points=[absVi(x,0,gamma_dot),absVi(x,1,gamma_dot),absVi(x,Sol(:)',gamma_dot)];
else
    Crit_points=[absVi(x,0,gamma_dot),absVi(x,1,gamma_dot)];
end

c=max(Crit_points)-Omega_max*r;
