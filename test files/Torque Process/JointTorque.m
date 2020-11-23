clear all
clc

load('Stance.mat');
load('Swing_Alex.mat');
% load('Swing_N1.mat')
% load('Swing_P2.mat')
% load('Newest.mat')

current_size = size(TorqueStance,1);
% num_of_samples = 50; %number of data points we wish to have per step

% TorqueSwing = interp1(linspace(0,1,num_of_samples),SwingTorque,linspace(0,1,current_size));
% TorqueSwing = interp1(linspace(0,1,num_of_samples),SwingTorqueN1,linspace(0,1,current_size));
% AnkleSwing = interp1(linspace(0,1,num_of_samples),SwingTorqueP2,linspace(0,1,current_size));
% TorqueSwing(:,2) = -TorqueSwing(:,2);
% TorqueSwing(:,3) = AnkleSwing(:,3);

T = {'Hip','Knee','Ankle'};

figure(1)
for i = 1:3
    subplot(3,1,i)
    plot(-TorqueStance(:,i))
    title(T{i})
end
set(gcf,'Position',[200 100 600 600])
% suptitle('Stance Torque from Emanuel')

figure(2)
for i = 1:3
    subplot(3,1,i)
    plot(TorqueSwing(:,i))
    title(T{i})
end
set(gcf,'Position',[200 100 600 600])

%% Smooth out data
TorqueALL = [TorqueSwing;-TorqueStance;TorqueSwing;-TorqueStance];

%number of blocks to average around
avgblocks = 3;
%converted to coefficients for the filtfilt function
coeffblocks = ones(1,avgblocks)/avgblocks;

%Smooth the dynamic data just like the kinematic data
for j=1:3
    for i=1:50
        TorqueALL(:,j) = filtfilt(coeffblocks,1,TorqueALL(:,j));
    end
end

TorqueALL =interp1(linspace(0,1,4*current_size),TorqueALL,linspace(0,1,4*num_of_samples));
Torque = TorqueALL(51:150,:);

figure(3)
for i = 1:3
    subplot(3,1,i)
    plot(Torque(:,i))
    title(T{i})
end