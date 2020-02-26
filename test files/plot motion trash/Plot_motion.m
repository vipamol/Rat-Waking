% Plot bodies and joints motion during walking
clear 
clc

load('motion.mat')

% Position of the joints
xh = pos_hip(:,1);
yh = pos_hip(:,2);
zh = pos_hip(:,3);
xk = pos_knee(:,1);
yk = pos_knee(:,2);
zk = pos_knee(:,3);
xa = pos_ankle(:,1);
ya = pos_ankle(:,2);
za = pos_ankle(:,3);

% Posion of the bodies (COM)
x1 = pos_femur(:,1);
y1 = pos_femur(:,2);
z1 = pos_femur(:,3);
x2 = pos_tibia(:,1);
y2 = pos_tibia(:,2);
z2 = pos_tibia(:,3);
x3 = pos_foot(:,1);
y3 = pos_foot(:,2);
z3 = pos_foot(:,3);

to_plot_2d = 0;
to_plot_3d = 1;


%% 2D plots
if to_plot_2d == 1
    figure(1) % X-Y motion
    for i = 1:length(dt) 
        plot([xh(i),xk(i)],[yh(i),yk(i)],'k-','Linewidth',2);
        hold on 
        plot([xk(i),xa(i)],[yk(i),ya(i)],'k-','Linewidth',2);
        hold on
        plot([xa(i),x3(i)],[ya(i),y3(i)],'k-','Linewidth',2);
        hold off
        
        title('X-Y motion(Side View)')
        axis([-0.04 0.04 -0.075 0])
        set(gca,'XDir','reverse');
        pbaspect([1 1 1])
        
        text(xh(i),yh(i),'o','color','r');
        text(x1(i),y1(i),'x','color','g')
        text(xk(i),yk(i),'o','color','r');
        text(x2(i),y2(i),'x','color','g');
        text(xa(i),ya(i),'o','color','r');
        text(x3(i),y3(i),'x','color','g');
        M(i) = getframe;
    end
%      movie(M,5,24)

    figure(2) % Z-x motion
    for i = 1:length(dt)
        plot([zh(i),zk(i)],[xh(i),xk(i)],'k-','Linewidth',2);
        hold on 
        plot([zk(i),za(i)],[xk(i),xa(i)],'k-','Linewidth',2);
        hold on 
        plot([za(i),z3(i)],[xa(i),x3(i)],'k-','Linewidth',2);
        hold off
        
        title('Z-X motion(Top View)')
        axis([-0.03 0 -0.04 0.04])
        pbaspect([1 1 1])
        
        text(zh(i),xh(i),'o','color','r');
        text(z1(i),x1(i),'x','color','g')
        text(zk(i),xk(i),'o','color','r');
        text(z2(i),x2(i),'x','color','g');
        text(za(i),xa(i),'o','color','r');
        text(z3(i),x3(i),'x','color','g');
        N(i) = getframe;
    end
%      movie(N,5,24)
    
    figure(3) %Y-Z motion
    for i = 1:length(dt)
        plot([zh(i),zk(i)],[yh(i),yk(i)],'k-','Linewidth',2);
        hold on 
        plot([zk(i),za(i)],[yk(i),ya(i)],'k-','Linewidth',2);
        hold on 
        plot([za(i),z3(i)],[ya(i),y3(i)],'k-','Linewidth',2);
        hold off
        
        title('Y-Z motion (Rear View)')
        axis([-0.03 0 -0.075 0])
        pbaspect([1 1 1])
        
        text(zh(i),yh(i),'o','color','r');
        text(z1(i),y1(i),'x','color','g')
        text(zk(i),yk(i),'o','color','r');
        text(z2(i),y2(i),'x','color','g');
        text(za(i),ya(i),'o','color','r');
        text(z3(i),y3(i),'x','color','g');
        O(i) = getframe;
    end
%      movie(O,5,24)
end


 %% 3D plots
    if to_plot_3d == 1
        figure(4) %3D motion
        for i = 1:length(dt)
            plot3([zh(i),zk(i)],[xh(i),xk(i)],[yh(i),yk(i)],'k-','Linewidth',2);
            hold on 
            plot3([zk(i),za(i)],[xk(i),xa(i)],[yk(i),ya(i)],'k-','Linewidth',2);
            hold on 
            plot3([za(i),z3(i)],[xa(i),x3(i)],[ya(i),y3(i)],'k-','Linewidth',2);
            hold off
            
            title('3d motion')
            axis([-0.03 0 -0.04 0.04 -0.075 0])
            pbaspect([1 1 1])
            grid on
           
            text(zh(i),xh(i),yh(i),'o','color','r');
            text(z1(i),x1(i),y1(i),'x','color','g')
            text(zk(i),xk(i),yk(i),'o','color','r');
            text(z2(i),x2(i),y2(i),'x','color','g');
            text(za(i),xa(i),ya(i),'o','color','r');
            text(z3(i),x3(i),y3(i),'x','color','g');
            L(i) = getframe;
        end
    end
%       movie(L,5,24)