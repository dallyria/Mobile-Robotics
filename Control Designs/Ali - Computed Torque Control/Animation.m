function Animation()
load Results4Anim
load Rob_param
l1=b; l2=b; Lo1=Lo; Lo2=Lo;
% d=0.17; %Robots diameter
% l12=0.12; %link btw. robot

figure
File_name=['TwoRobotsDynamics','.gif'];

% des_traj=plot(x_ref,y_ref,'-.b');

x_c = q(:,1);
y_c = q(:,2);
theta1 = q(:,3);
theta2 = q(:,4);

l=1;
for i=1:length(x_c)
    x_1=x_c(i);
    y_1=y_c(i);
    theta_1=theta1(i);
Poi1=[x_1+2*l1*cos(theta_1+pi/4)/sqrt(2),y_1+2*l1*sin(theta_1+pi/4)/sqrt(2);...
    x_1+2*l1*cos(theta_1+pi/4+pi/2)/sqrt(2),y_1+2*l1*sin(theta_1+pi/4+pi/2)/sqrt(2);...
    x_1+2*l1*cos(theta_1+pi/4+pi)/sqrt(2),y_1+2*l1*sin(theta_1+pi/4+pi)/sqrt(2);...
    x_1+2*l1*cos(theta_1+pi/4+3*pi/2)/sqrt(2),y_1+2*l1*sin(theta_1+pi/4+3*pi/2)/sqrt(2);
    x_1+2*l1*cos(theta_1+pi/4)/sqrt(2),y_1+2*l1*sin(theta_1+pi/4)/sqrt(2)];

    x_2=x_c2(i);
    y_2=y_c2(i);
    theta_2=theta2(i);
    
Poi2=[x_2+2*l2*cos(theta_2+pi/4)/sqrt(2),y_2+2*l2*sin(theta_2+pi/4)/sqrt(2);...
    x_2+2*l2*cos(theta_2+pi/4+pi/2)/sqrt(2),y_2+2*l2*sin(theta_2+pi/4+pi/2)/sqrt(2);...
    x_2+2*l2*cos(theta_2+pi/4+pi)/sqrt(2),y_2+2*l2*sin(theta_2+pi/4+pi)/sqrt(2);...
    x_2+2*l2*cos(theta_2+pi/4+3*pi/2)/sqrt(2),y_2+2*l2*sin(theta_2+pi/4+3*pi/2)/sqrt(2);
    x_2+2*l2*cos(theta_2+pi/4)/sqrt(2),y_2+2*l2*sin(theta_2+pi/4)/sqrt(2)];

robot1=plot(Poi1(:,1),Poi1(:,2),'r');
hold on
%   wheels 1st robot:
left1=plot([x_1+(m+mo)*2*l1*cos(theta_1+pi/4)/sqrt(2),x_1+0.6*2*l1*cos(theta_1+pi/4+pi/2)/sqrt(2)],...
    [y_1+(m+mo)*2*l1*sin(theta_1+pi/4)/sqrt(2),y_1+(m+mo)*2*l1*sin(theta_1+pi/4+pi/2)/sqrt(2)],'b');

right1=plot([x_1+(m+mo)*2*l1*cos(theta_1+pi/4+pi)/sqrt(2),x_1+0.6*2*l1*cos(theta_1+pi/4+3*pi/2)/sqrt(2)],...
    [y_1+(m+mo)*2*l1*sin(theta_1+pi/4+pi)/sqrt(2),y_1+(m+mo)*2*l1*sin(theta_1+pi/4+3*pi/2)/sqrt(2)],'m');

%   wheels 2nd robot:
left2=plot([x_2+(m+mo)*2*l2*cos(theta_2+pi/4)/sqrt(2),x_2+(m+mo)*2*l2*cos(theta_2+pi/4+pi/2)/sqrt(2)],...
    [y_2+(m+mo)*2*l2*sin(theta_2+pi/4)/sqrt(2),y_2+(m+mo)*2*l2*sin(theta_2+pi/4+pi/2)/sqrt(2)],'b');

right2=plot([x_2+(m+mo)*2*l2*cos(theta_2+pi/4+pi)/sqrt(2),x_2+(m+mo)*2*l2*cos(theta_2+pi/4+3*pi/2)/sqrt(2)],...
    [y_2+(m+mo)*2*l2*sin(theta_2+pi/4+pi)/sqrt(2),y_2+(m+mo)*2*l2*sin(theta_2+pi/4+3*pi/2)/sqrt(2)],'m');

% %   Docking Link:
Dock1=[x_1+0.5*2*l1*cos(theta_1+pi),y_1+0.5*2*l1*sin(theta_1+pi);...
    x_1+(0.5*2*l1+Lo1)*cos(theta_1+pi),y_1+(0.5*2*l1+Lo1)*sin(theta_1+pi)];

Dock2=[x_2+2*l2*cos(theta_2)/2,y_2+2*l2*sin(theta_2)/2;...
    x_2+(0.5*2*l2+Lo2)*cos(theta_2),y_2+(0.5*2*l2+Lo2)*sin(theta_2)];


dock1=plot(Dock1(:,1),Dock1(:,2),'r');
dock2=plot(Dock2(:,1),Dock2(:,2),'k');
plot(x_2+(0.5*2*l2+Lo2)*cos(theta_2),y_2+(0.5*2*l2+Lo2)*sin(theta_2),'ko','MarkerSize',2)
robot2=plot(Poi2(:,1),Poi2(:,2),'k');
legend([robot1,robot2],{'Leader','Follower'})

%% ---- Axis
%axis([x_1-1 x_1+1 y_1-1 y_1+1]) % Moving axis

axis([-2 2 -2 2]) % fix axis
%%


axis equal
grid on
xlabel x(m)
ylabel y(m)
title('dynamic motion of robots')

text(-.25+x_1,.5+y_1,['Time =',num2str(T(i),3),'(sec)'])
 f = getframe;
 im=frame2im(f);
 [imind,cm] = rgb2ind(im,256);



    if i == 1;
          imwrite(imind,cm,File_name,'gif', 'Loopcount',inf);
    else
       imwrite(imind,cm,File_name,'gif','WriteMode','append'); 
    end
    hold off

       l=l+1;
end


