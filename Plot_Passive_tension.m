function Plot_Passive_tension()
% Plot muscle passive tension during walking
load('Passive_tension.mat')

Plot_all = 0; % if Plot_all = 1 plot all muscle together, otherwise plot seperately.

num = length(muscle_name);
for i = 1:num
    muscle_name{i} = strrep(muscle_name{i},'LH_','');
end


if Plot_all ==1
    for i = 1:num
        plot(Tension{i},'Linewidth',1.5)
        Legend{i} = muscle_name{i};
        hold on 
    end
    title('Muscle Passive Tension during walking')
    set(gcf,'Position',[50 100 700 500])
    legend(Legend,'FontSize',8,'Position',[0.45 0.6 0 0])    %In case you want show the legend
else
    for i = 1:num
        subplot(2,num/2,i)
        plot(Tension{i},'Linewidth',1.5)
        title(muscle_name{i},'FontSize',9)
    end
    set(gcf,'Position',[200 100 1200 500])
    suptitle('Muscle Passive Tension during walking')
end

end