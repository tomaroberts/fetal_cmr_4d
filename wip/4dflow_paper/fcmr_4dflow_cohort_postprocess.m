%% Analysis: 4D Flow Paper Post-processing
% script to keep track of and to maintain consistent processing between 
% fcmr subjects used in paper.

%% studyDir
studyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper';


%% fcmr Numbers for paper:

% 4 stacks in recon
% fcmr4stacks = [206 254];

% >5 stacks in recon
fcmr5stacks = [189 191 194 197 201 202 214];


%% fcmr189
fNum = 189; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

fcmr_4dflow_postprocessing( studyDir, fNum );
fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr' );


%% fcmr191
fNum = 191; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

fcmr_4dflow_postprocessing( studyDir, fNum );
fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr' );


%% fcmr194
fNum = 194; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

% default
bloodpoolMask = 'mask_blood_pool_cleaned';
fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr', 'bloodpoolMask', bloodpoolMask );

%%% MASKS
% All:
velMasks = {'mask_aorta','mask_LV','mask_RV','mask_LOT','mask_ROT','mask_LA_RA','mask_IVC_SVC'};
% % Right-sided
% velMasks = {'mask_RV','mask_ROT','mask_LA_RA','mask_IVC_SVC'};
% % Left-sided
% velMasks = {'mask_aorta','mask_LV','mask_LOT'};
% % Left-sided w/out LOT
% velMasks = {'mask_aorta','mask_LV'};

bloodpoolMask = 'mask_blood_pool_cleaned';
fcmr_4dflow_postprocessing( studyDir, fNum, ...
    'bloodpoolMask', bloodpoolMask, ...
    'velMasks', velMasks);
fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr', ...
    'bloodpoolMask', bloodpoolMask, ...
    'velMasks', velMasks);


%% fcmr197
fNum = 197; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

fcmr_4dflow_postprocessing( studyDir, fNum );
fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr' );


%% fcmr201
fNum = 201; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

fcmr_4dflow_postprocessing( studyDir, fNum );
fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr' );

% % newer SVRTK recon with -transformations option
% % accidentally used -slice_transformations option in prior recon
% % - have now become ../vel_vol and ../vel_vol_4d
% fcmr_4dflow_postprocessing( studyDir, fNum, 'velDir', 'vel_vol_trans' );
% fcmr_4dflow_postprocessing( studyDir, fNum, 'velDir', 'vel_vol_trans', 'fileExt', 'polyCorr' );


%% fcmr202
studyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr202';
fNum = 202; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

fcmr_4dflow_postprocessing( studyDir, fNum );
fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr' );

%%% MASKS
% All
velMasks = {'mask_aorta','mask_LV','mask_RV','mask_LOT','mask_ROT','mask_LA_RA','mask_IVC_SVC','mask_PA_DA'};

fcmr_4dflow_postprocessing( studyDir, fNum, ...
    'velMasks', velMasks);
fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr', ...
    'velMasks', velMasks);

% % newer SVRTK recon with -transformations option
% % accidentally used -slice_transformations option in prior recon
% % - have now become ../vel_vol and ../vel_vol_4d
% fcmr_4dflow_postprocessing( studyDir, fNum, 'velDir', 'vel_vol_trans', ...
%     'velMasks', velMasks);
% fcmr_4dflow_postprocessing( studyDir, fNum, 'velDir', 'vel_vol_trans', 'fileExt', 'polyCorr', ...
%     'velMasks', velMasks);

%%% v2 - removes some erroneous vectors in RPA and lower DAo
%%% MASKS
% All
velMasks = {'mask_aorta_v2','mask_LV','mask_RV','mask_LOT','mask_ROT','mask_LA_RA','mask_IVC_SVC','mask_PA_DA_v2'};

% fcmr_4dflow_postprocessing( studyDir, ...
%     'velMasks', velMasks);
fcmr_4dflow_postprocessing( studyDir, 'fileExt', 'polyCorr', ...
    'velMasks', velMasks);



%% fcmr214
fNum = 214; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

fcmr_4dflow_postprocessing( studyDir, fNum );
fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr' );


%% c_fcmr246 - ToF
studyDir = 'I:\fcmr_4d_chloe_recons\c_fcmr246';
fNum = 246; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

% default
% bloodpoolMask = 'mask_blood_pool_cleaned';
% fcmr_4dflow_postprocessing( studyDir, fNum, 'fileExt', 'polyCorr', 'bloodpoolMask', bloodpoolMask );

%%% MASKS
% All:
% nb: basically folded aorta into LVOT in this instance, rather than having
% mask_aorta
velMasks = {'mask_LVOT','mask_RVOT','mask_RV','mask_LV','mask_LA_RA','mask_IVC_SVC'};

% Ao / PA only
% velMasks = {'Ao_Milou','PA_Milou'};


bloodpoolMask = 'mask_blood_pool_cleaned';
% fcmr_4dflow_postprocessing( studyDir, fNum, ...
%     'bloodpoolMask', bloodpoolMask, ...
%     'velMasks', velMasks);
fcmr_4dflow_postprocessing( studyDir, 'fileExt', 'polyCorr', ...
    'bloodpoolMask', bloodpoolMask, ...
    'velMasks', velMasks);


%% c_fcmr240_tom - CAT
studyDir = 'I:\fcmr_4d_chloe_recons\c_fcmr240_tom';

fNum = 240; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

% default
% bloodpoolMask = 'mask_blood_pool_cleaned';
% fcmr_4dflow_postprocessing( studyDir, 'fileExt', 'polyCorr', 'bloodpoolMask', bloodpoolMask );

% polyCorr
bloodpoolMask = 'mask_blood_pool_cleaned';
velMasks      = bloodpoolMask;
fcmr_4dflow_postprocessing( studyDir, 'useVelDriftCorr', true, 'fileExt', 'polyCorr', 'bloodpoolMask', bloodpoolMask, 'velMasks', velMasks );


%%% MASKS
% All:
velMasks = {'mask_CAT','mask_Ao','mask_DAo','mask_DA','mask_IVC_SVC','mask_LA_RA','mask_LV','mask_RV','mask_PA'};


bloodpoolMask = 'mask_blood_pool_cleaned';
% fcmr_4dflow_postprocessing( studyDir, fNum, ...
%     'bloodpoolMask', bloodpoolMask, ...
%     'velMasks', velMasks);
fcmr_4dflow_postprocessing( studyDir, 'fileExt', 'polyCorr', ...
    'bloodpoolMask', bloodpoolMask, ...
    'velMasks', velMasks);


%% c_fcmr329
studyDir = 'I:\King's College London\Perinatal CMR - Documents\General\Data\fcmr_4d_clinical\c_fcmr329';

fNum = 329; disp(['Running fcmr ' num2str(fNum) ' ...' ]);

% default
% bloodpoolMask = 'mask_blood_pool_cleaned';
% fcmr_4dflow_postprocessing( studyDir, 'fileExt', 'polyCorr', 'bloodpoolMask', bloodpoolMask );

% polyCorr
bloodpoolMask = 'mask_blood_pool_cleaned';
velMasks      = bloodpoolMask;
fcmr_4dflow_postprocessing( studyDir, 'useVelDriftCorr', true, 'fileExt', 'polyCorr', 'bloodpoolMask', bloodpoolMask, 'velMasks', velMasks );


%%% MASKS
% All:
velMasks = {'AAO','DAO','LSVC','PA','PA_window','SVC'};


bloodpoolMask = 'mask_blood_pool_cleaned';
% fcmr_4dflow_postprocessing( studyDir, fNum, ...
%     'bloodpoolMask', bloodpoolMask, ...
%     'velMasks', velMasks);
fcmr_4dflow_postprocessing( studyDir, 'fileExt', 'polyCorr', ...
    'bloodpoolMask', bloodpoolMask, ...
    'velMasks', velMasks);