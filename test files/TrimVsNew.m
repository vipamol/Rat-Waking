clear all
clc
load('JointKinematics2.mat')
Theta = Theta-pi/2;
Theta_old = Theta;
current_size_old = size(Theta_old,1);

load('JointKinematics.mat')
Theta = Theta-pi/2;
Theta_new = Theta;
current_size_new = size(Theta_new,1);

% Theta_old = [Theta_old;Theta_old;Theta_old];
% Theta_new = [Theta_new;Theta_new;Theta_new];

num_of_samples = 100;
Theta_old = interp1(linspace(0,1,current_size_old),Theta_old,linspace(0,1,num_of_samples));
Theta_new = interp1(linspace(0,1,current_size_new),Theta_new,linspace(0,1,num_of_samples));

figure(1)
plot(Theta_new(:,1),'b:','Linewidth',2)
hold on
plot(Theta_old(:,1),'r:','Linewidth',2)
hold off
title('Hip Angles (rads)')

figure(2)
plot(Theta_new(:,2),'b:','Linewidth',2)
hold on
plot(Theta_old(:,2),'r:','Linewidth',2)
hold off
title('Knee Angles (rads)')

figure(3)
plot(Theta_new(:,3),'b:','Linewidth',2)
hold on
plot(Theta_old(:,3),'r:','Linewidth',2)
hold off
title('Ankle Angles (rads)')