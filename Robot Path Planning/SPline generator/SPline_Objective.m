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


function F= SPline_Objective(x)

global Points

tP=Points(:,3);
xP=Points(:,1);
yP=Points(:,2);

for i=1:size(Points,1)
   r=Ri(x,tP(i));
   E(i,1)=r(1)-xP(i);
   E(i,2)=r(2)-yP(i);
end
F=sum(sum(E.^2,2));