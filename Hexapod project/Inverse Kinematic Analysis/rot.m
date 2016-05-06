%##################################################################
%
%
% This program is about to find a trajectory of Transitional Maryland 
% University Manipulator.
% This is an Assignment of Advanced Robotics Course -
% Sharif University of Technology - International Campus
% Student ID : Shahriari 89251844

function RT=rot(phi)

RT=[cos(phi),sin(phi),0;-sin(phi),cos(phi),0;0,0,1];