L = 35:0.1:51.5;
Lw = 51.5-35;
Lr1 = 51.5;
Lr2 = 35 + Lw/2;

for i = 1:166
    
Y1(i) = (1-(L(i)-Lr1)^2/(Lw^2))*100;
Y2(i) = (1-(L(i)-Lr2)^2/(Lw^2))*100;
end

subplot(1,2,1)
plot(L,Y1)
ylim([-5 105])
xlim([34.5 52])
xlabel('(mm)')
ylabel('% Isometric Tension  Used')
title('Y = (1-(L-Lrest)^2/Lwidth^2)*100')

subplot(1,2,2)
plot(L,Y2)
ylim([70 102])
xlim([34.5 52])
xlabel('(mm)')
ylabel('% Isometric Tension  Used')

title('Y = (1-(L-Lrest)^2/Lwidth^2)*100')