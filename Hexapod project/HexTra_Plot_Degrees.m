% This program Plots the trajectory of HEXAPOD ROBOT
% ADVANCED ROBOTIC PROJECT - Supervisor Dr Osguie
% Shahriari 89251844


%% getting the data of solved trajectory
clear t_L1 t_L2 t_L3 t_L4 t_L5 t_L6
TA=180*zero22pi(t_array,'radian')/pi;

for L=1:size(t_array,2)/6
t_L1(:,L)=TA(1:3,6*L-5);
t_L2(:,L)=TA(1:3,6*L-4);
t_L3(:,L)=TA(1:3,6*L-3);
t_L4(:,L)=TA(1:3,6*L-2);
t_L5(:,L)=TA(1:3,6*L-1);
t_L6(:,L)=TA(1:3,6*L);
end

figure(2)
subplot(6,3,1)
plot(time_array,(t_L1(1,:)))
ylabel '\Theta_C_1'
title 'Coxa'
subplot(6,3,2)
plot(time_array,t_L1(2,:))
ylabel '\Theta_f_1'
title 'Femur'
subplot(6,3,3)
plot(time_array,t_L1(3,:))
ylabel '\Theta_t_1'
title 'Tibia'

subplot(6,3,4)
plot(time_array,(t_L2(1,:)))
ylabel '\Theta_C_2'
subplot(6,3,5)
plot(time_array,t_L2(2,:))
ylabel '\Theta_f_2'
subplot(6,3,6)
plot(time_array,t_L2(3,:))
ylabel '\Theta_t_2'

subplot(6,3,7)
plot(time_array,t_L3(1,:))
ylabel '\Theta_C_3'
subplot(6,3,8)
plot(time_array,t_L3(2,:))
ylabel '\Theta_f_3'
subplot(6,3,9)
plot(time_array,t_L3(3,:))
ylabel '\Theta_t_3'

subplot(6,3,10)
plot(time_array,t_L4(1,:))
ylabel '\Theta_C_4'
subplot(6,3,11)
plot(time_array,t_L4(2,:))
ylabel '\Theta_f_4'
subplot(6,3,12)
plot(time_array,t_L4(3,:))
ylabel '\Theta_t_4'


subplot(6,3,13)
plot(time_array,t_L5(1,:))
ylabel '\Theta_C_5'
subplot(6,3,14)
plot(time_array,t_L5(2,:))
ylabel '\Theta_f_5'
subplot(6,3,15)
plot(time_array,t_L5(3,:))
ylabel '\Theta_t_5'

subplot(6,3,16)
plot(time_array,t_L6(1,:))
ylabel '\Theta_C_6'
xlabel 'Time (sec)'

subplot(6,3,17)
plot(time_array,t_L6(2,:))
ylabel '\Theta_f_6'
xlabel 'Time (sec)'

subplot(6,3,18)
plot(time_array,t_L6(3,:))
ylabel '\Theta_t_6'
xlabel 'Time (sec)'


