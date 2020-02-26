load('Johnson2011MuscleParameters.mat')
Musc_par = cell(9,size(params,2));

for i = 1:size(params,2)
    Musc_par{1,i} = params{1,i};   %Names
    Musc_par{2,i} = params{5,i};   %BFA
    Musc_par{3,i} = params{6,i};   %BFP
    Musc_par{4,i} = params{42,i};  %VA
    Musc_par{5,i} = params{22,i};  %GA
    Musc_par{6,i} = params{40,i};  %TA
    Musc_par{7,i} = params{38,i};  %SO
    Musc_par{8,i} = params{20,i};  %IP
    Musc_par{9,i} = params{34,i};  %RF
end

save('BiarMuscleParameters.mat','Musc_par')