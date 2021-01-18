%% Wrapper script for vessel analysis in 4D flow paper
% used to generate flow curves fed into Excel/Graphpad

studyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper';


%% 7 used in paper
fcmr5stacks = [189 191 194 197 201 202 214];

for ff = 1:numel(fcmr5stacks)
    fcmr_4dflow_vessel_analysis( studyDir, fcmr5stacks(ff), false )
    fcmr_4dflow_vessel_analysis( studyDir, fcmr5stacks(ff), true )
end


%% Other 2 reconstructions
fcmr4stacks = [206 254];
for ff = 1:numel(fcmr4stacks)
    fcmr_4dflow_vessel_analysis( studyDir, fcmr4stacks(ff), false )
    fcmr_4dflow_vessel_analysis( studyDir, fcmr4stacks(ff), true )
end