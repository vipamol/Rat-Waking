%Smooth Stride smooths out the mean animal kinematic data for the hind legs
%a bit and then writes it, the derivative, and the 2nd derivative to a
%Signal file for use in with the SimMechanics simulation of the hind leg
%which calculates the torques to produce those movements. Need to run
%CompareAnimalAndSim first.

clear thisstridec thistime

%How many data points to average
avgblocks = 20;
coeffblocks = ones(1,avgblocks)/avgblocks;

%Convert the mean angle trace to radians
thisstride = AllAnimalMean*pi/180;

%Grab just the swing part
thistime = Time(Time>.5)';
n=size(thisstride,1);
m = size(thistime,1);

%Crop the angle data to the swing part
thisstridec = thisstride(n-m+1:end,:);

%Make the first time at 0
thistime = thistime-thistime(1);

figure
plot(thistime,thisstridec)
hold on

%Make the vector to length of t
t=1000;
thistimec = linspace(thistime(1),thistime(end),t);
thisstridec = interp1(thistime,thisstridec,thistimec);

for i=1:150
    thisstridec = filtfilt(coeffblocks,1,thisstridec);
end

% figure
plot(thistimec,thisstridec)

%Differentiate
dydt = diff(thisstridec)*t;

figure
plot(thistimec(1:end-1),dydt)

%Differentiate again
dydt2 = diff(dydt)*t;

figure
plot(thistimec(1:end-2),dydt2)

%Write angle and derivatives to files, make sure they are all the same
%length
% xlswrite('C:\Users\Kaiyu\Desktop\Rat Kinematics and Dynamics\Kinematic and Dynamic Processing\Signal.xlsx',[thistimec',thisstridec],1)
% xlswrite('C:\Users\Kaiyu\Desktop\Rat Kinematics and Dynamics\Kinematic and Dynamic Processing\Signal.xlsx',[thistimec(1:end-1)',dydt],2)
% xlswrite('C:\Users\Kaiyu\Desktop\Rat Kinematics and Dynamics\Kinematic and Dynamic Processing\Signal.xlsx',[thistimec(1:end-2)',dydt2],3)