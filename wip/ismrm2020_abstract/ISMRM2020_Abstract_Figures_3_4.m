%% Analysis: 4D Flow Paper Post-processing
% script to keep track of and to maintain consistent processing between 
% fcmr subjects

% This script is for the HRH reconstructions
% Made for abstract submission ISMRM 2020


%% studyDir
studyDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fcmr194 - Normal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fcmr194 --- standard k-t

reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr194';
cd(reconDir);

fNum = 194; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

% default
bloodpoolMask = 'mask_blood_pool_cleaned';

%%% MASKS
% All:
velMasks = {'mask_aorta','mask_LV','mask_RV','mask_LOT','mask_ROT','mask_LA_RA','mask_IVC_SVC'};

fcmr_4dflow_postprocessing( reconDir, 'fileExt', 'polyCorr', ...
    'bloodpoolMask', bloodpoolMask, ...
    'velMasks', velMasks);


%% fcmr194_hrh_fullRecon

reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr194_hrh_fullRecon';
cd(reconDir);

fNum = 194; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

% default
bloodpoolMask = 'mask_blood_pool_slightlyCleaned';
% fcmr_4dflow_postprocessing( reconDir, 'fileExt', 'polyCorr', 'bloodpoolMask', bloodpoolMask );

%%% MASKS
% All:
velMasks = {'mask_aorta','mask_LV','mask_RV','mask_LOT','mask_ROT','mask_LA_RA','mask_IVC_SVC'};

fcmr_4dflow_postprocessing( reconDir, 'fileExt', 'polyCorr', ...
    'bloodpoolMask', bloodpoolMask, ...
    'velMasks', velMasks);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fcmr202 - RAA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% fcmr202 - standard k-t

reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr202';
cd(reconDir);

fNum = 202; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

%%% MASKS
% All
bloodpoolMask = 'mask_blood_pool';
velMasks = {'mask_aorta','mask_LV','mask_RV','mask_LOT','mask_ROT','mask_LA_RA','mask_IVC_SVC','mask_PA_DA'};

fcmr_4dflow_postprocessing( reconDir, 'fileExt', 'polyCorr', ...
    'bloodpoolMask', bloodpoolMask, ...
    'velMasks', velMasks);


%% fcmr202 - hrh k-t

reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr202_hrh_fullRecon';
cd(reconDir);

fNum = 202; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

%%% MASKS
% All
bloodpoolMask = 'mask_blood_pool_cleaned';
velMasks = {'mask_aorta','mask_LV','mask_RV','mask_LOT','mask_ROT','mask_LA_RA','mask_IVC_SVC','mask_PA_DA'};

fcmr_4dflow_postprocessing( reconDir, 'fileExt', 'polyCorr', ...
    'bloodpoolMask', bloodpoolMask, ...
    'velMasks', velMasks);

