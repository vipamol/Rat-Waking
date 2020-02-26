%Unfortunately we do not have data that is connected for stance and swing
%kinematics and dynamics, so we must get a little creative in how we put it
%all together. We calculate swing torques from Manuela's kinematic data on
%limb compliance input into a SimMechanic model of the rat. We then mash
%the stance and swing data together and smooth it out.

% function output = ConnectData()

%First read in all the Rat data from excel and mat files
[ATime, AFootContact, AFrontLeft, AFrontRight, ABackLeft, ABackRight] = Process_Ratte_Kinematics;
load('mean_kinematics_stand.mat');
load('mean_torques.mat')

%This cycles through the legs (currently set to just the hind legs) finds
%when the steps occured, breaks them up, puts stance as the first 50% and
%swing as the second 50% and stores this data away
[Stride, StrideContact, NewTime] = Animal_Stride_Normalizing(0, ATime, AFootContact, AFrontLeft, AFrontRight, ABackLeft, ABackRight);
% load('Animal_Stride_Normalizing.mat')

%The part of the code used here is to get the mean stride kinematic angles
%and plots it
Ani_RMS;

%This smooths out the swing data and finds the derivatives. This data is
%used for simulating in the SimMechanic Model
SmoothStride

%This finds the hip, knee, and ankle angles based on the input kinematic
%data of femur, fibula, and foot and assumes a fixed pelvis angle of 30
%degrees for stance data
Hip_theta = pi/6+femur_theta_0MW;
Knee_theta = pi-fibula_theta_0MW+femur_theta_0MW;
Ankle_theta = pi-fibula_theta_0MW+foot_theta_0MW;
KinematicsStance = [Hip_theta' Knee_theta' Ankle_theta'];
save('KinematicsStance.mat','KinematicsStance')

%thisstridec is a value output by the SmoothStride function
KinematicsSwing = thisstridec;

%Torque data is stored into a vector and converted from the normalized form
%into the Newton-meters and resaved
Torque0 = [Hip_0MW', Knee_0MW', Ankle_0MW'];
TorqueinNm = Torque0*9.81*.30112;
save('TorqueinNm.mat','TorqueinNm');

% TorqueStancestruct = load('C:\Users\Hunt\Desktop\PhD\Optimizers\Torque Calculation\TorqueinNm.mat');
% KinematicsStancestruct = load('C:\Users\Hunt\Desktop\PhD\Optimizers\Torque Calculation\KinematicsStance.mat');

% TorqueStance = TorqueStancestruct.TorqueinNm;
% KinematicsStance = KinematicsStancestruct.KinematicsStance;

%Crop off bad data (NaNs) and flip torque direction of Knee (it is in
%opposite direction of hip and ankle)
KinematicsStance = KinematicsStance(1:end-5,:);
TorqueStance = TorqueinNm(1:end-5,:);
TorqueStance(:,2) = -TorqueStance(:,2);

%Record lengths for use in normalizing data between 0 and 1
m = length(TorqueStance);
n = length(KinematicsSwing);

%Make kinematic swing data as long as stance data and put it all together
KinematicsSwing = interp1(1:n,KinematicsSwing,linspace(1,n,m));
KinematicsAll = [KinematicsStance;KinematicsSwing];
t=size(KinematicsAll,1);

%This is for plotting the data on top of the previous kinematic data
%commented out here and done after smoothing and final data is made
for i=1:3
    figure(i)
    hold on
    plot(0:1/(2*m):1-1/(2*m),KinematicsAll(:,i)*180/pi,'--')
end

close all

%% Smoothing

%number of blocks to average around
avgblocks = 10;
%converted to coefficients for the filtfilt function
coeffblocks = ones(1,avgblocks)/avgblocks;

%Pad the data so begining of stance and end of swing match
KinematicsAll2 = [KinematicsAll;KinematicsAll;KinematicsAll];

%For each joint, do the averaging of 10 data points many (originally set to
%50) times
for j=1:3
    for i=1:50
        KinematicsAll2(:,j) = filtfilt(coeffblocks,1,KinematicsAll2(:,j));
    end
end

%Crop off the padded data - Kinematic Data complete
KinematicsAll = KinematicsAll2(t+1:2*t,:);

%Plot the data on the other kinematic data figures
for i=1:3
    figure(i)
    hold on
    plot(0:1/(2*m):1-1/(2*m),KinematicsAll(:,i)*180/pi,'o:')
end

%Differentiate the full data for use in the 'optimizer' code
Athisstridec = KinematicsAll2;
Adydt = diff(Athisstridec)*t/0.28;
Adydt2 = diff(Adydt)*t/0.28;

%Crop the data to the same length
Theta = KinematicsAll;
Theta_dot = Adydt(t+1:2*t,:);
Theta_doubledot = Adydt2(t+1:2*t,:);

%Save the data
AllKinematicInfo = [.28*[0:1/t:1-1/t]',Theta,Theta_dot,Theta_doubledot]';
save('AllSignal.mat','AllKinematicInfo')
save('JointKinematics','Theta','Theta_dot','Theta_doubledot')
keyboard
%% Use new kinematic data to find torques on joints during swing

%Get just the swing data part and differentiate it twice for simulation
thisstridec = KinematicsAll(t/2:end,:);
dydt = diff(thisstridec)*t/0.28;
dydt2 = diff(dydt)*t/0.28;

%Put all the data into one vector for saving and using the SimMechanics
%model
AllSwingKinematicInfo = [.28*[0:1/t:.5-2/t]',thisstridec(1:end-2,:),dydt(1:end-1,:),dydt2]';
save('Signal.mat','AllSwingKinematicInfo')

%Simulate the system with the developed swing data
sim('HindLegSwing.slx')

%Store the timeseries vectors generated from the simulation into vectors
Ankle_T = AnkleTorque.data;
Knee_T = KneeTorque.data;
Hip_T = HipTorque.data;

%%

%Put all the joints together and make it the same length as the kinematic
%data for swing
TorqueSwing = [Hip_T Knee_T Ankle_T];
% TorqueSwing = TorqueSwing(1:end-1,:);
l = length(TorqueSwing);
TorqueSwing = interp1(1:l,TorqueSwing,linspace(1,l,m));

%Pad the data for smoothing
TorqueAll2 = [TorqueSwing;-TorqueStance;TorqueSwing;-TorqueStance];

%Plot for fun
h = figure;
plotvec = linspace(0,1,length(TorqueAll2));
plot(plotvec,1e3*TorqueAll2(:,:),'linewidth',2)
title('Joint Torque')
xlabel('Percent Stride')
ylabel('Torque (mN-m)')
set(h,'Position',[500,500,700,250])
grid on


%number of blocks to average around
avgblocks = 3;
%converted to coefficients for the filtfilt function
coeffblocks = ones(1,avgblocks)/avgblocks;

%Smooth the dynamic data just like the kinematic data
for j=1:3
    for i=1:50
        TorqueAll2(:,j) = filtfilt(coeffblocks,1,TorqueAll2(:,j));
    end
end

%Crop out the excess that was added for padding and save
TorqueAll = TorqueAll2(m+1:3*m-2,:);
sizetorque = size(TorqueAll);
sizetheta = size(Theta);
save('AllTraining2_6_20.mat','TorqueAll','Theta','Theta_dot')
save('Kinematics and Dynamics.mat','TorqueAll','Theta','Theta_dot','Theta_doubledot')



%Plot for fun
h = figure;
plotvec = linspace(0,1,length(TorqueAll));
plot(plotvec,1e3*TorqueAll(:,:),'linewidth',2)
title('Joint Torque')
xlabel('Percent Stride')
ylabel('Torque (mN-m)')
set(h,'Position',[500,500,700,250])
grid on

%Plot for fun
h = figure;
plotvec = linspace(0,1,length(KinematicsAll));
plot(plotvec,KinematicsAll(:,:),'linewidth',2)
title('Joint Angle')
xlabel('Percent Stride')
ylabel('Angle (Radians)')
set(h,'Position',[500,500,700,250])
grid on