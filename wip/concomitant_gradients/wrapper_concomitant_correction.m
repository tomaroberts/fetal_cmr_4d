
reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr194';
cd(reconDir);
cd data

sID = 16;
sliceToDisplay = 5;
scaleFactor = 1;

% reslice_nii( ['s' num2str(sID) '_rlt_ph.nii.gz'], ['s' num2str(sID) '_rlt_ph-RESLICE.nii.gz'] );

fcmr_concomitant_correction( reconDir, sID, ...
    'sliceToDisplay', sliceToDisplay, 'scaleFactor', scaleFactor );