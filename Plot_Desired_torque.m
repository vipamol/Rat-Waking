function Plot_Desired_torque()

load('Kinematics and Dynamics.mat')
TorqueAll = interp1(linspace(0,1,588),TorqueAll,linspace(0,1,100));

% load('EOM.mat')
% Desired = TorqueAll-EOM;

T = {'Hip';'Knee';'Ankle'};

for i = 1:3
    subplot(1,3,i)
    plot(TorqueAll(:,i),'b','Linewidth',1.5)
    title(T{i});
end
suptitle('Desired Joint Torques')

keyboard
end