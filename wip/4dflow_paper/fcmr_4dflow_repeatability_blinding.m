%% Randomisation of fcmr_4d_flow data for David and Milou
% - analysis performed at request of Nat Comms reviewers
% - intention: intra-/inter-user stats and reliability/confidence analysis

studyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper';
cd(studyDir);

analysisDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\#4dflow_paper_analysis_info\flow_analysis_repeatability';
cine_volFilename = 'cine_vol-RESLICE.nii.gz';

paperIDs = [189 191 194 197 201 202 214 240 246 308];
numCases = numel(paperIDs);


%% Copy, paste cine_vol files and rename based on randomisation

randID = randperm(numCases*2);
randID = reshape(randID,numCases,2);

for ff = 1:numel(paperIDs)
      
    if paperIDs(ff) == 240
        cd('I:\fcmr_4d_chloe_recons\c_fcmr240_tom');
    elseif paperIDs(ff) == 246
        cd('I:\fcmr_4d_chloe_recons\c_fcmr246');        
    elseif paperIDs(ff) == 308
        cd('I:\fcmr_4d_chloe_recons\c_fcmr308');
    else
        cd( fullfile( studyDir, ['fcmr' num2str(paperIDs(ff)) ] ) );
    end
    
    cd cine_vol
    
    % Copy twice for repeatability
    copyfile( cine_volFilename, [analysisDir '\' num2str(randID(ff,1)) '_cine_vol.nii.gz'] );
    copyfile( cine_volFilename, [analysisDir '\' num2str(randID(ff,2)) '_cine_vol.nii.gz'] );
    
end


%% Save IDs for unblinding
fcmrIDs = [paperIDs' randID];

cd(analysisDir);
save('fcmrIDs_for_unblinding','fcmrIDs','-v7');