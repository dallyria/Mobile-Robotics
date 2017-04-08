% This program shows the trajectory animation of HEXAPOD ROBOT
% ADVANCED ROBOTIC PROJECT - Supervisor Dr Osguie
% Shahriari 89251844


for K=1:size(Joints_GF_Array,2)/4 %Trajectory


%% getting the data
Joints_GF=Joints_GF_Array(1:18,4*K-3:4*K);% Make an array of all points in each configuration for every 4 coloumns.
%% Plot Hexapod===========================================================
figure(1)
Pthyln=-20:10:TotStp*Vy*delta_t/steps+20;
Pthxln=zeros(size(Pthyln));
Pthzln=zeros(size(Pthyln));
plot3(Pthxln,Pthyln,Pthzln,'-->b','MarkerEdgeColor','r','MarkerSize',7)

hold on
plot3(20,-20,0,-20,0,20,'w')
for i=1:6
    h=plot3(Joints_GF(3*i-2,:),Joints_GF(3*i-1,:),Joints_GF(3*i,:),'-ro','Linewidth',2,...
        'MarkerEdgeColor','k',...
                'MarkerFaceColor','c',...
                'MarkerSize',6);
end

for i=1:6
h=plot3([Joints_GF(1,1),Joints_GF(4,1),Joints_GF(7,1),...
    Joints_GF(10,1),Joints_GF(13,1),Joints_GF(16,1),Joints_GF(1,1)],...
    [Joints_GF(2,1),Joints_GF(5,1),Joints_GF(8,1),...
    Joints_GF(11,1),Joints_GF(14,1),Joints_GF(17,1),Joints_GF(2,1)],...
    [Joints_GF(3,1),Joints_GF(6,1),Joints_GF(9,1),...
    Joints_GF(12,1),Joints_GF(15,1),Joints_GF(18,1),Joints_GF(3,1)],...
    'm','Linewidth',3);
    
end
hold off
axis equal
axis manual
grid on
BJ=[[Joints_GF(1,1),Joints_GF(4,1),Joints_GF(7,1),...
    Joints_GF(10,1),Joints_GF(13,1),Joints_GF(16,1),Joints_GF(1,1)];...
    [Joints_GF(2,1),Joints_GF(5,1),Joints_GF(8,1),...
    Joints_GF(11,1),Joints_GF(14,1),Joints_GF(17,1),Joints_GF(2,1)];...
    [Joints_GF(3,1),Joints_GF(6,1),Joints_GF(9,1),...
    Joints_GF(12,1),Joints_GF(15,1),Joints_GF(18,1),Joints_GF(3,1)]];% Body joints in ground frame

text(BJ(1,1),BJ(2,1),BJ(3,1)+1,'#1')
text(BJ(1,2),BJ(2,2),BJ(3,2)+1,'#2')
text(BJ(1,3),BJ(2,3),BJ(3,3)+1,'#3')
text(BJ(1,4),BJ(2,4),BJ(3,4)+1,'#4')
text(BJ(1,5),BJ(2,5),BJ(3,5)+1,'#5')
text(BJ(1,6),BJ(2,6),BJ(3,6)+1,'#6')

xlabel ' X(cm) '
ylabel ' Y(cm) '
zlabel ' Z(cm) '   
title 'Trajectory of Six legged Walking Robot'
F(K)=getframe;
% G(K)=getframe(gcf); % Snapshut the plot
end

%% Showing the animation of the TRAJECTORY

movie(F,2,10)
% movie2avi(G, 'hex_mov.avi','fps',100,'quality',70);