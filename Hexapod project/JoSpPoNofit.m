p11=polyfit(time_array,t_L1(1,:),7);
t_L1pn(1,:)=polyval(p11,time_array);
plot(time_array,t_L1(1,:),'--r')
hold on
plot(time_array,t_L1pn(1,:),'k')
hold off

p12=polyfit(time_array,t_L1(2,:),7);
t_L1pn(2,:)=polyval(p12,time_array);
plot(time_array,t_L1(2,:),'--r')
hold on
plot(time_array,t_L1pn(2,:),'k')
hold off