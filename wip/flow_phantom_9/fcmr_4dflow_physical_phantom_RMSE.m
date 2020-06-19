
cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9');
load('voxelwise_velocity_values_fromGraphpad','AP','FH','RL','MAG');

% In paper:
% AP = y
% RL = x
% FH = z
% MAG = MAG


% AP = [];
% FH = [];
% RL = [];
% MAG = [];

% copy columns from Graphpad
% These values generated in C:\Users\tr17\Documents\CODE_Projects\fetal_cmr_4dflow\wip\flow_phantom_9\fcmr_4dflow_fp9_comparison.m
% But couldn't get it working for some reason... think I adjusted code and
% need to re-examine

MEAN.AP = ( mean( AP(:) ) );
MEAN.FH = ( mean( FH(:) ) );
MEAN.RL = ( mean( RL(:) ) );
MEAN.MAG = ( mean( MAG(:) ) );

MSE.AP = ( immse( AP(:,1), AP(:,2) ) );
MSE.FH = ( immse( FH(:,1), FH(:,2) ) );
MSE.RL = ( immse( RL(:,1), RL(:,2) ) );
MSE.MAG = ( immse( MAG(:,1), MAG(:,2) ) );

RMSE.AP = sqrt( immse( AP(:,1), AP(:,2) ) );
RMSE.FH = sqrt( immse( FH(:,1), FH(:,2) ) );
RMSE.RL = sqrt( immse( RL(:,1), RL(:,2) ) );
RMSE.MAG = sqrt( immse( MAG(:,1), MAG(:,2) ) );

%%% nMSE

% denominator = range
nMSE.AP = 100 * (MSE.AP / (max(AP(:)) - min(AP(:))));
nMSE.FH = 100 * (MSE.FH / (max(FH(:)) - min(FH(:))));
nMSE.RL = 100 * (MSE.RL / (max(RL(:)) - min(RL(:))));
nMSE.MAG = 100 * (MSE.MAG / (max(MAG(:)) - min(MAG(:))));

% denominator = 99% percentile
nMSEprctile.AP = 100 * (MSE.AP / ( prctile(AP, 99, 'all') ) );
nMSEprctile.FH = 100 * (MSE.FH / ( prctile(FH, 99, 'all') ) );
nMSEprctile.RL = 100 * (MSE.RL / ( prctile(RL, 99, 'all') ) );
nMSEprctile.MAG = 100 * (MSE.MAG / ( prctile(MAG, 99, 'all') ) );

% denominator = mean
nMSEmean.AP = 100 * (MSE.AP / mean(abs(AP(:))) );
nMSEmean.FH = 100 * (MSE.FH / mean(abs(FH(:))) );
nMSEmean.RL = 100 * (MSE.RL / mean(abs(RL(:))) );
nMSEmean.MAG = 100 * (MSE.MAG / mean(MAG(:)) );


%%% nRMSE

% denominator = range
nRMSE.AP = 100 * (RMSE.AP / (max(AP(:)) - min(AP(:))));
nRMSE.FH = 100 * (RMSE.FH / (max(FH(:)) - min(FH(:))));
nRMSE.RL = 100 * (RMSE.RL / (max(RL(:)) - min(RL(:))));
nRMSE.MAG = 100 * (RMSE.MAG / (max(MAG(:)) - min(MAG(:))));

% denominator = 99% percentile
nRMSEprctile.AP = 100 * (RMSE.AP / ( prctile(AP, 99, 'all') ) );
nRMSEprctile.FH = 100 * (RMSE.FH / ( prctile(FH, 99, 'all') ) );
nRMSEprctile.RL = 100 * (RMSE.RL / ( prctile(RL, 99, 'all') ) );
nRMSEprctile.MAG = 100 * (RMSE.MAG / ( prctile(MAG, 99, 'all') ) );

% denominator = mean
nRMSEmean.AP = 100 * (RMSE.AP / mean(abs(AP(:))) );
nRMSEmean.FH = 100 * (RMSE.FH / mean(abs(FH(:))) );
nRMSEmean.RL = 100 * (RMSE.RL / mean(abs(RL(:))) );
nRMSEmean.MAG = 100 * (RMSE.MAG / mean(MAG(:)) );






% % think this doesn't work because range -ve to +ve
% nRMSEmean.AP = 100 * (RMSE.AP / mean(abs(AP(:))) );
% nRMSEmean.FH = 100 * (RMSE.FH / mean(abs(FH(:))) );
% nRMSEmean.RL = 100 * (RMSE.RL / mean(abs(RL(:))) );
% nRMSEmean.MAG = 100 * (RMSE.MAG / mean(MAG(:)) );
% 
% nRMSEstd.AP = 100 * (RMSE.AP / std(AP(:)) );
% nRMSEstd.FH = 100 * (RMSE.FH / std(FH(:)) );
% nRMSEstd.RL = 100 * (RMSE.RL / std(RL(:)) );
% nRMSEstd.MAG = 100 * (RMSE.MAG / std(MAG(:)) );


