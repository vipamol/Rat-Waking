clear all
clc

load('motion.mat')
load('attach.mat')

% 1-BFA & 7-IP(Hip); 2_BFP & 8-RF(Hip&Knee);  5-TA & 6-SO(Ankle); 
% 3-VA & 4-GA (Knee, GA^ankle)
muscle1 = [musc_pos{1,1} musc_pos{1,2} musc_pos{1,3} musc_pos{1,4}];
muscle2 = [musc_pos{2,1} musc_pos{2,2} musc_pos{2,3}];
muscle3 = [musc_pos{3,1} musc_pos{3,2} musc_pos{3,3} musc_pos{3,4}];
muscle4 = [musc_pos{4,1} musc_pos{4,2} musc_pos{4,3}];
muscle5 = [musc_pos{5,1} musc_pos{5,2} musc_pos{5,3} musc_pos{5,4}];
muscle6 = [musc_pos{6,1} musc_pos{6,2} musc_pos{6,3}];
muscle7 = [musc_pos{7,1} musc_pos{7,2} musc_pos{7,3}];
muscle8 = [musc_pos{8,1} musc_pos{8,2} musc_pos{8,3} musc_pos{8,4}];

% Wanna see hip muscles with motion, Plot 1,2,7,8;
% Wanna see knee muscles with motion, Plot 2,3,4,8;
% Wanna see ankle muscles with motion, Plot 4,5,6;

for i = 1:length(dt)
    % plot the joint and body
            plot3([pos_hip(i,3),pos_knee(i,3)],[pos_hip(i,1),pos_knee(i,1)],[pos_hip(i,2),pos_knee(i,2)],'k-','Linewidth',3);
            hold on 
            plot3([pos_knee(i,3),pos_ankle(i,3)],[pos_knee(i,1),pos_ankle(i,1)],[pos_knee(i,2),pos_ankle(i,2)],'k-','Linewidth',3);
            hold on 
            plot3([pos_ankle(i,3),pos_foot(i,3)],[pos_ankle(i,1),pos_foot(i,1)],[pos_ankle(i,2),pos_foot(i,2)],'k-','Linewidth',3);
            hold on
            
   % plot muscles   1     
            plot3([muscle1(i,6),muscle1(i,3)],[muscle1(i,4),muscle1(i,1)],[muscle1(i,5),muscle1(i,2)],'r--','Linewidth',1);
            hold on 
            plot3([muscle1(i,9),muscle1(i,6)],[muscle1(i,7),muscle1(i,4)],[muscle1(i,8),muscle1(i,5)],'r--','Linewidth',1);
            hold on 
            plot3([muscle1(i,12),muscle1(i,9)],[muscle1(i,10),muscle1(i,7)],[muscle1(i,11),muscle1(i,8)],'r--','Linewidth',1);
            hold on 
   % plot muscles   2     
            plot3([muscle2(i,6),muscle2(i,3)],[muscle2(i,4),muscle2(i,1)],[muscle2(i,5),muscle2(i,2)],'b--','Linewidth',1);
            hold on 
            plot3([muscle2(i,9),muscle2(i,6)],[muscle2(i,7),muscle2(i,4)],[muscle2(i,8),muscle2(i,5)],'b--','Linewidth',1);
            hold on 
   % plot muscles   3     
            plot3([muscle3(i,6),muscle3(i,3)],[muscle3(i,4),muscle3(i,1)],[muscle3(i,5),muscle3(i,2)],'g--','Linewidth',1);
            hold on 
            plot3([muscle3(i,9),muscle3(i,6)],[muscle3(i,7),muscle3(i,4)],[muscle3(i,8),muscle3(i,5)],'g--','Linewidth',1);
            hold on 
            plot3([muscle3(i,12),muscle3(i,9)],[muscle3(i,10),muscle3(i,7)],[muscle3(i,11),muscle3(i,8)],'g--','Linewidth',1);
            hold on 
   % plot muscles   4     
            plot3([muscle4(i,6),muscle4(i,3)],[muscle4(i,4),muscle4(i,1)],[muscle4(i,5),muscle4(i,2)],'c--','Linewidth',1);
            hold on 
            plot3([muscle4(i,9),muscle4(i,6)],[muscle4(i,7),muscle4(i,4)],[muscle4(i,8),muscle4(i,5)],'c--','Linewidth',1);
            hold on 
   % plot muscles   5     
            plot3([muscle5(i,6),muscle5(i,3)],[muscle5(i,4),muscle5(i,1)],[muscle5(i,5),muscle5(i,2)],'m--','Linewidth',1);
            hold on 
            plot3([muscle5(i,9),muscle5(i,6)],[muscle5(i,7),muscle5(i,4)],[muscle5(i,8),muscle5(i,5)],'m--','Linewidth',1);
            hold on 
            plot3([muscle5(i,12),muscle5(i,9)],[muscle5(i,10),muscle5(i,7)],[muscle5(i,11),muscle5(i,8)],'m--','Linewidth',1);
            hold on 
   % plot muscles   6     
            plot3([muscle6(i,6),muscle6(i,3)],[muscle6(i,4),muscle6(i,1)],[muscle6(i,5),muscle6(i,2)],'m--','Linewidth',1);
            hold on 
            plot3([muscle6(i,9),muscle6(i,6)],[muscle6(i,7),muscle6(i,4)],[muscle6(i,8),muscle6(i,5)],'m--','Linewidth',1);
            hold on 
   %  plot muscles   7     
            plot3([muscle7(i,6),muscle7(i,3)],[muscle7(i,4),muscle7(i,1)],[muscle7(i,5),muscle7(i,2)],'r--','Linewidth',1);
            hold on 
            plot3([muscle7(i,9),muscle7(i,6)],[muscle7(i,7),muscle7(i,4)],[muscle7(i,8),muscle7(i,5)],'r--','Linewidth',1);
            hold on 
   % plot muscles   8     
            plot3([muscle8(i,6),muscle8(i,3)],[muscle8(i,4),muscle8(i,1)],[muscle8(i,5),muscle8(i,2)],'b--','Linewidth',1);
            hold on 
            plot3([muscle8(i,9),muscle8(i,6)],[muscle8(i,7),muscle8(i,4)],[muscle8(i,8),muscle8(i,5)],'b--','Linewidth',1);
            hold on 
            plot3([muscle8(i,12),muscle8(i,9)],[muscle8(i,10),muscle8(i,7)],[muscle8(i,11),muscle8(i,8)],'b--','Linewidth',1);
            hold off 
           
   
            title('3d motion')
            axis equal
            axis([-0.03 0 -0.05 0.05 -0.07 0])
            pbaspect([1 1 1])
          
%             set(gcf,'Position',[400 100 700 600])  
%              view([180 0])  % Any View


             %%% 3 Views
%             view([-90 0])   % Side View
%             view([-180 0])  % Front view
%             view([0 90])    % Top View
            grid on
           
            
%             % joints and bodies
%             text(pos_hip(i,3),pos_hip(i,1),pos_hip(i,2),'o','color','r');
%             text(pos_femur(i,3),pos_femur(i,1),pos_femur(i,2),'x','color','g')
%             text(pos_knee(i,3),pos_knee(i,1),pos_knee(i,2),'o','color','r');
%             text(pos_tibia(i,3),pos_tibia(i,1),pos_tibia(i,2),'x','color','g');
%             text(pos_ankle(i,3),pos_ankle(i,1),pos_ankle(i,2),'o','color','r');
%             text(pos_foot(i,3),pos_foot(i,1),pos_foot(i,2),'x','color','g');
            M(i) = getframe;
end

% movie(M,5,24)