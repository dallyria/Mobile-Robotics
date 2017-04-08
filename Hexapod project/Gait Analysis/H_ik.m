%%This Program calculates the inverse kinematic of Hexapod's one limb
%%Shahriari
%%ADVANCED ROBOTIC COURSE PROJECT
%---------------------------------------------------------------------

function  T= H_ik(X)
x=X(1);
y=X(2);
z=X(3);
X_o=[0,0,0];
%hexapod_limb_inverse_kinematic
l0=2.5;%cm
l1=7.4;%cm
l2=11.4;%cm

t1=atan2(y,x)
c=((x-l0*cos(t1))^2+(y-l0*sin(t1))^2+z^2)^.5;
B=acos((-l2^2+l1^2+c^2)/(2*l1*c));
t2=asin(-z/c)-B
C1=pi/2-B;
h=l1*sin(B);
C2=acos(h/l2);
t3=pi-C1-C2
% T=[1,t1;2,t2;3,t3]';
% save q_input.mat T

x1=X_o(1);y1=X_o(2);z1=X_o(3);
x2=l0*cos(t1);y2=l0*sin(t1);z2=0;
x3=(l0+l1*cos(t2))*cos(t1);
y3=(l0+l1*cos(t2))*sin(t1);
z3=-l1*sin(t2);
x4=(l0+l2*cos(t3+t2)+l1*cos(t2))*cos(t1);
y4=(l0+l2*cos(t3+t2)+l1*cos(t2))*sin(t1);
z4=-l1*sin(t2)-l2*sin(t3+t2);
x4
y4 
z4
x_a=[x1 x2 x3 x4];
y_a=[y1 y2 y3 y4];
z_a=[z1 z2 z3 z4];
for i=1:3
plot3([x_a(i) x_a(i+1)],[y_a(i) y_a(i+1)],[z_a(i) z_a(i+1)],'k','linewidth',2)
hold on
end
xlabel 'Position x(cm)'
ylabel 'Position z(cm)'
zlabel 'Position y(cm)'
grid on
axis equal
hold off

