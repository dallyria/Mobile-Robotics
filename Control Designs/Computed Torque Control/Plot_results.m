close all
warning('off','all')
load MultiRoboSysParameters

sim('Ali_model')
q=Q_out_sim;
figure
% Leader's Trajectory
L=plot(q(:,1),q(:,2));
hold on
plot(q(1,1),q(1,2),'*');
% Follower's Trajectory
x_c2 = q(:,1)-h*cos(q(:,3))-h*cos(q(:,4));
y_c2 = q(:,2)-h*sin(q(:,3))-h*sin(q(:,4));
F=plot(x_c2,y_c2,'r');
plot(x_c2(1),y_c2(1),'r*');

% Desired Trajectory
clear q_desired

qd=desired_trajectory';
D=plot(qd(1,:),qd(2,:),'k');
legend([L F D],{'Leader trajectory','Follower Trajectory','Desired trajectory'})
hold off
grid on
xlabel x(m);
ylabel y(m);
%axis([-.5 2.5 0 1.4])
%axis equal

% plot errors:
figure
simTime=linspace(0,tout(end),length(error));
subplot 411
plot(simTime,error(:,1))
title('Error signals')
ylabel('x')
xlabel('time')
grid on

subplot 412
plot(simTime,error(:,2))
ylabel('y')
xlabel('time')
grid on
subplot 413
plot(simTime,error(:,3))
ylabel('\theta_1')
xlabel('time')
grid on
subplot 414
plot(simTime,error(:,4))
ylabel('\theta_2')
xlabel('time')
grid on
% 
% clear tau
% tau=tau_oufromSIM1;
% figure()
% plot(tau(:,1));
% hold
% plot(tau(:,2),'r');
% plot(tau(:,3),'g');
% plot(tau(:,4),'k');
% legend('Leader Left Wheel','Leader Right Wheel','Follower Left Wheel','Follower Right Wheel');
T=simTime;
save Results4Anim T q x_c2 y_c2;    
% Animation;