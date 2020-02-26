function Plot_muscle()
% Plot muscle length and velocity during walking
load('muscle.mat')

Plot_all = 0; % if Plot_all = 1 plot all muscle together, otherwise plot seperately.

num = length(muscle_name);
for i = 1:num
    muscle_name{i} = strrep(muscle_name{i},'LH_','');
end
Legend = cell(muscle_name);

%% plot muscle length
figure(1)
if Plot_all ==1
    for i = 1:num
        plot(muscle_length{i},'Linewidth',1.5)
        Legend{i} = muscle_name{i};
        hold on 
    end
    title('Muscle Length during walking')
    set(gcf,'Position',[50 100 700 500])
    legend(Legend,'FontSize',8,'Position',[0.5 0.4 0 0])    %In case you want show the legend
else
    for i = 1:num
        subplot(2,num/2,i)
        plot(muscle_length{i},'Linewidth',1.5)
        title(muscle_name{i},'FontSize',9)
    end
    set(gcf,'Position',[200 100 1200 500])
    suptitle('Muscle Length during walking')
end

%% plot muscle velocity
figure(2)
if Plot_all ==1
    for i = 1:num
        plot(muscle_velocity{i},'Linewidth',1.5)
        Legend{i} = muscle_name{i};
        hold on 
    end
    title('Muscle Velocity during walking')
    set(gcf,'Position',[800 100 700 500])
    legend(Legend,'FontSize',8,'Position',[0.77 0.25 0 0])    %In case you want show the legend
else
    for i = 1:num
        subplot(2,num/2,i)
        plot(muscle_velocity{i},'Linewidth',1.5)
        title(muscle_name{i},'FontSize',9)
    end
    set(gcf,'Position',[200 100 1200 500])
    suptitle('Muscle velocity during walking')
end

end