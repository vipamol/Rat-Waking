%Find RMS for animal and simulation data
close all

Time = 0:.001:1;

for k=1:3

clear AnimalStride AnimalMean AnimalStd SimStride SimMean SimStd X Y
length(Stride)
for i=1:length(Stride)
        AnimalStride(:,i) = interp1(NewTime{i,3},Stride{i,3}(:,k),Time);
end

for j=1:size(AnimalStride,1)
    clear temp
    temp=AnimalStride(j,:);
    AnimalMean(j) = mean(temp(~isnan(temp)),2);
    AnimalStd(j) = std(temp(~isnan(temp)));

end

figure(k)

hold on
X=[Time,fliplr(Time)];
Y=[AnimalMean+AnimalStd,fliplr(AnimalMean-AnimalStd)];
fill(X,Y,[.5 .5 1])
plot(Time,AnimalMean,'--b','Linewidth',2)

AllAnimalMean(:,k) = AnimalMean';

end


%Begin Front Legs

for k=1:4

clear AnimalStride AnimalMean AnimalStd SimStride SimMean SimStd  X Y

for i=1:length(Stride)
    if isempty(Stride{i,2})==0
        AnimalStride(:,i) = interp1(NewTime{i,2},Stride{i,2}(:,k),Time);
    else
       AnimalStride(:,i) = Time*NaN; 
    end
end

for j=1:size(AnimalStride,1)
    clear temp
    temp=AnimalStride(j,:);
    AnimalMean(j) = mean(temp(~isnan(temp)),2);
    AnimalStd(j) = std(temp(~isnan(temp)));
end

figure
hold on
X=[Time,fliplr(Time)];
Y=[AnimalMean+AnimalStd,fliplr(AnimalMean-AnimalStd)];
fill(X,Y,[.5 .5 1])
plot(Time,AnimalMean,'--b','Linewidth',2)

end

