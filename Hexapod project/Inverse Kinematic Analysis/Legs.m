function LEGS_inv = Legs(t)
global Ltip_JF l0 l1 l2 lt lc lf 

t1=t(1,1:6)';
t2=t(2,1:6)';
t3=t(3,1:6)';

for i=1:6    
    x=Ltip_JF(1,i);
    y=Ltip_JF(2,i);
    z=Ltip_JF(3,i);
       
        LEGS_inv(3*i-2)=t1(i)-atan2(y,x);
        c=((x-l0*cos(t1(i)))^2+(y-l0*sin(t1(i)))^2+z^2)^.5;
        B=acos((-l2^2+l1^2+c^2)/(2*l1*c));
        LEGS_inv(3*i-1)=t2(i)-(asin(-z/c)-B);
        C1=pi/2-B;
        h=l1*sin(B);
        C2=acos(h/l2);
        LEGS_inv(3*i)=t3(i)-(pi-C1-C2);
        
%         LEGS_inv(3*i-2)=cos(t1(i))*(lc+lf*cos(t2(i))+lt*cos(t2(i)+t3(i)))-x;
%         LEGS_inv(3*i-1)=sin(t1(i))*(lc+lf*cos(t2(i))+lt*cos(t2(i)+t3(i)))-y;
%         LEGS_inv(3*i)=-lf*sin(t2(i))-lt*sin(t3(i)+t2(i))-z;
%         
end

