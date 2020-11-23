clear all
clc

load('Muscle_ten.mat');

%number of blocks to average around
avgblocks = 20;
%converted to coefficients for the filtfilt function
coeffblocks = ones(1,avgblocks)/avgblocks;

%Smooth the dynamic data just like the kinematic data
for j=1:8
    for i=1:50
        TensionALL(:,j) = filtfilt(coeffblocks,1,Muscle_tension(:,j));
    end
end

for i = 1:100
    if TensionALL(i,1) < 0
        TensionALL(i,7) = -TensionALL(i,1);
        TensionALL(i,1) = 0;
    end
    
    if TensionALL(i,2) < 0
        TensionALL(i,8) = -TensionALL(i,2);
        TensionALL(i,2) = 0;
    end
    
    if TensionALL(i,4) < 0
        TensionALL(i,3) = -TensionALL(i,4);
        TensionALL(i,4) = 0;
    end
    
    if TensionALL(i,5) < 0
        TensionALL(i,6) = -TensionALL(i,5);
        TensionALL(i,5) = 0;
    end
end

T = {'BFA','BFP','VA','GA','TA','SO','IP','RF'};

for i = 1:8
    subplot(2,4,i)
    plot(TensionALL(:,i),'Linewidth',1.5)
    title(T{i})
end
set(gcf,'Position',[200 100 1200 500])
suptitle('Tension(Ib feedback) during Walking')