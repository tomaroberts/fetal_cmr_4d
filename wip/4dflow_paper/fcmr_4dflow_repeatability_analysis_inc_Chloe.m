%% fcmr_4dflow_repeatability_analysis
% - March 2020 - NatComms reviewers requested reliability analysis.
% - Blinded data created using: E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\#4dflow_paper_analysis_info\flow_analysis_repeatability\fcmr_4dflow_repeatability_blinding.m
% - David/Milou = expert ROI drawers
%
%


%% Directories
studyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper';
blindingDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\#4dflow_paper_analysis_info\flow_analysis_repeatability';


%% Unblind
cd(blindingDir);
load('fcmrIDs_for_unblinding.mat');

paperIDs     = fcmrIDs([1,2,3,4,5,6,7],1); 
numCases     = numel(paperIDs);
fetalWeights = [2.11, 1.82, 0.88, 1.48, 1.82, 1.99, 1.01]; % kg, taken from 4D_Flow_Paper_Fetal_volumes_from_Milou.xls

trial1IDs = fcmrIDs([1,2,3,4,5,6,7],2);
trial2IDs = fcmrIDs([1,2,3,4,5,6,7],3); 

roiDirPath = @( ii, roiDirName ) fullfile( studyDir, ['fcmr' num2str(paperIDs(ii))], roiDirName );


%% Copy ROIs to studyDir

% % Milou
% for ii = 1:numCases
%     cd( fullfile( blindingDir, 'analysis_milou', num2str(trial1IDs(ii)) ) );
%     copyfile( '*.nii.gz' , roiDirPath( ii, 'roi_milou_trial1' ) );   
%     cd( fullfile( roiDirPath( ii, 'roi_milou_trial1' ) ) ); delete('*cine_vol.nii.gz'); % ugly, but easiest
%     
%     cd( fullfile( blindingDir, 'analysis_milou', num2str(trial2IDs(ii)) ) );
%     copyfile( '*.nii.gz' , roiDirPath( ii, 'roi_milou_trial2' ) );
%     cd( fullfile( roiDirPath( ii, 'roi_milou_trial2' ) ) ); delete('*cine_vol.nii.gz'); % ugly, but easiest
% end
% 
% % David
% for ii = 1:numCases
%     cd( fullfile( blindingDir, 'analysis_david', num2str(trial1IDs(ii)) ) );
%     copyfile( '*.nii.gz' , roiDirPath( ii, 'roi_david_trial1' ) );   
%         cd( fullfile( roiDirPath( ii, 'roi_david_trial1' ) ) ); delete('*cine_vol.nii.gz'); % ugly, but easiest
%     
%     cd( fullfile( blindingDir, 'analysis_david', num2str(trial2IDs(ii)) ) );
%     copyfile( '*.nii.gz' , roiDirPath( ii, 'roi_david_trial2' ) );
%     cd( fullfile( roiDirPath( ii, 'roi_david_trial2' ) ) ); delete('*cine_vol.nii.gz'); % ugly, but easiest
% end


%%% Vessel analysis
%% Milou
for ff = 1:numCases
    
    roiDir = 'roi_milou_trial1';
    disp(['Running ' num2str(paperIDs(ff)) ' Trial 1 ...' ]);
    T1 = fcmr_4dflow_vessel_analysis( studyDir, paperIDs(ff), true, roiDir ); % true = polyCorr
    T1.FMAG = T1.FMAG(:, 1:ceil(size(T1.FMAG,2)/2) );                         % ugly 1:... = ignore STD columns in FMAG
    
    roiDir = 'roi_milou_trial2';
    disp(['Running ' num2str(paperIDs(ff)) ' Trial 2 ...' ]);
    T2 = fcmr_4dflow_vessel_analysis( studyDir, paperIDs(ff), true, roiDir ); % true = polyCorr
    T2.FMAG = T2.FMAG(:, 1:ceil(size(T2.FMAG,2)/2) );                         % ugly 1:... = ignore STD columns in FMAG
    
    % collate flow value tables
    v = strcat('fcmr',num2str(paperIDs(ff)));
    MILOU.(v).T1 = T1.FMAG;
    MILOU.(v).T2 = T2.FMAG;
    
    % convert flow to ml/min/kg
    MILOU.(v).T1{:,2:end} = ( T1.FMAG{:,2:end} .* 60 ) ./ fetalWeights(ff);   % ml/min/kg
    MILOU.(v).T2{:,2:end} = ( T2.FMAG{:,2:end} .* 60 ) ./ fetalWeights(ff);   % ml/min/kg
    
    % calculate mean/std flow
    MILOU.(v).T1mean = mean( MILOU.(v).T1{:,2:end}, 1 );
    MILOU.(v).T1std  =  std( MILOU.(v).T1{:,2:end}, 1 );
    
    MILOU.(v).T2mean = mean( MILOU.(v).T2{:,2:end}, 1 );
    MILOU.(v).T2std  =  std( MILOU.(v).T2{:,2:end}, 1 );
    
    % add fetalWeights to stucture for saving
    MILOU.(v).fetalWeight = fetalWeights(ff);

end

save( fullfile( blindingDir, 'milou_flow_analysis.mat' ), 'MILOU' );
close all;
disp('Finished MILOU analysis ...');

cd(blindingDir);


%% David
for ff = 1:numCases
    
    roiDir = 'roi_david_trial1';
    disp(['Running ' num2str(paperIDs(ff)) ' Trial 1 ...' ]);
    T1 = fcmr_4dflow_vessel_analysis( studyDir, paperIDs(ff), true, roiDir ); % true = polyCorr
    T1.FMAG = T1.FMAG(:, 1:ceil(size(T1.FMAG,2)/2) );                         % ugly 1:... = ignore STD columns in FMAG
    
    roiDir = 'roi_david_trial2';
    disp(['Running ' num2str(paperIDs(ff)) ' Trial 2 ...' ]);
    T2 = fcmr_4dflow_vessel_analysis( studyDir, paperIDs(ff), true, roiDir ); % true = polyCorr
    T2.FMAG = T2.FMAG(:, 1:ceil(size(T2.FMAG,2)/2) );                         % ugly 1:... = ignore STD columns in FMAG
    
    % collate flow value tables
    v = strcat('fcmr',num2str(paperIDs(ff)));
    DAVID.(v).T1 = T1.FMAG;
    DAVID.(v).T2 = T2.FMAG;
    
    % convert flow to ml/min/kg
    DAVID.(v).T1{:,2:end} = ( T1.FMAG{:,2:end} .* 60 ) ./ fetalWeights(ff);   % ml/min/kg
    DAVID.(v).T2{:,2:end} = ( T2.FMAG{:,2:end} .* 60 ) ./ fetalWeights(ff);   % ml/min/kg
    
    % calculate mean/std flow
    DAVID.(v).T1mean = mean( DAVID.(v).T1{:,2:end}, 1 );
    DAVID.(v).T1std  =  std( DAVID.(v).T1{:,2:end}, 1 );
    
    DAVID.(v).T2mean = mean( DAVID.(v).T2{:,2:end}, 1 );
    DAVID.(v).T2std  =  std( DAVID.(v).T2{:,2:end}, 1 );
    
    % add fetalWeights to stucture for saving
    DAVID.(v).fetalWeight = fetalWeights(ff);

end

save( fullfile( blindingDir, 'david_flow_analysis.mat' ), 'DAVID' );
close all;
disp('Finished DAVID analysis ...');

cd(blindingDir);


%% Chloe
% nb: 1 trial only so far
for ff = 1:numCases
    
    roiDir = 'roi_chloe_trial1';
    disp(['Running ' num2str(paperIDs(ff)) ' Trial 1 ...' ]);
    T1 = fcmr_4dflow_vessel_analysis( studyDir, paperIDs(ff), true, roiDir ); % true = polyCorr
    T1.FMAG = T1.FMAG(:, 1:ceil(size(T1.FMAG,2)/2) );                         % ugly 1:... = ignore STD columns in FMAG
    
    % collate flow value tables
    v = strcat('fcmr',num2str(paperIDs(ff)));
    CHLOE.(v).T1 = T1.FMAG;
    
    % convert flow to ml/min/kg
    CHLOE.(v).T1{:,2:end} = ( T1.FMAG{:,2:end} .* 60 ) ./ fetalWeights(ff);   % ml/min/kg
    
    % calculate mean/std flow
    CHLOE.(v).T1mean = mean( CHLOE.(v).T1{:,2:end}, 1 );
    CHLOE.(v).T1std  =  std( CHLOE.(v).T1{:,2:end}, 1 );
    
    % add fetalWeights to stucture for saving
    CHLOE.(v).fetalWeight = fetalWeights(ff);

end

save( fullfile( blindingDir, 'chloe_flow_analysis.mat' ), 'CHLOE' );
close all;
disp('Finished CHLOE analysis ...');

cd(blindingDir);


%% Collate Mean Flows through Time
cd(blindingDir);

load( fullfile( blindingDir, 'milou_flow_analysis.mat' ), 'MILOU' );
load( fullfile( blindingDir, 'david_flow_analysis.mat' ), 'DAVID' );
load( fullfile( blindingDir, 'chloe_flow_analysis.mat' ), 'CHLOE' );


%% Mean Flows

%%% IMPORTANT: order here is different to Excel sheet!
%%% Copied these values into:
%%% fcmr_4DFlow_ROI_reliability_analysis_FINAL.xlsx
vesselNames = MILOU.fcmr189.T1.Properties.VariableNames(2:6);

% Milou
MILOU.meanFlow.T1 = nan(numCases,5);
MILOU.meanFlow.T2 = nan(numCases,5);
MILOU.stdFlow.T1  = nan(numCases,5);
MILOU.stdFlow.T2  = nan(numCases,5);

for ff = 1:numCases

    v = strcat('fcmr',num2str(paperIDs(ff)));
    
    if strcmp(v,'fcmr194') % no DA
        MILOU.meanFlow.T1(ff,[1,3:5]) = MILOU.(v).T1mean;
        MILOU.meanFlow.T2(ff,[1,3:5]) = MILOU.(v).T2mean;
        MILOU.stdFlow.T1(ff,[1,3:5]) = MILOU.(v).T1std;
        MILOU.stdFlow.T2(ff,[1,3:5]) = MILOU.(v).T2std;
    else
        MILOU.meanFlow.T1(ff,:) = MILOU.(v).T1mean;
        MILOU.meanFlow.T2(ff,:) = MILOU.(v).T2mean;
        MILOU.stdFlow.T1(ff,:) = MILOU.(v).T1std;
        MILOU.stdFlow.T2(ff,:) = MILOU.(v).T2std;
    end
end


% David
DAVID.meanFlow.T1 = nan(numCases,5);
DAVID.meanFlow.T2 = nan(numCases,5);
DAVID.stdFlow.T1  = nan(numCases,5);
DAVID.stdFlow.T2  = nan(numCases,5);

for ff = 1:numCases

    v = strcat('fcmr',num2str(paperIDs(ff)));
    
    if strcmp(v,'fcmr194') % no DA
        DAVID.meanFlow.T1(ff,[1,3:5]) = DAVID.(v).T1mean;
        DAVID.meanFlow.T2(ff,[1,3:5]) = DAVID.(v).T2mean;
        DAVID.stdFlow.T1(ff,[1,3:5]) = DAVID.(v).T1std;
        DAVID.stdFlow.T2(ff,[1,3:5]) = DAVID.(v).T2std;
    else
        DAVID.meanFlow.T1(ff,:) = DAVID.(v).T1mean;
        DAVID.meanFlow.T2(ff,:) = DAVID.(v).T2mean;
        DAVID.stdFlow.T1(ff,:) = DAVID.(v).T1std;
        DAVID.stdFlow.T2(ff,:) = DAVID.(v).T2std;
    end
end


% Chloe
CHLOE.meanFlow.T1 = nan(numCases,5);
CHLOE.stdFlow.T1  = nan(numCases,5);

for ff = 1:numCases

    v = strcat('fcmr',num2str(paperIDs(ff)));
    
    if strcmp(v,'fcmr194') % no DA
        CHLOE.meanFlow.T1(ff,[1,3:5]) = CHLOE.(v).T1mean;
        CHLOE.stdFlow.T1(ff,[1,3:5]) = CHLOE.(v).T1std;
    else
        CHLOE.meanFlow.T1(ff,:) = CHLOE.(v).T1mean;
        CHLOE.stdFlow.T1(ff,:) = CHLOE.(v).T1std;
    end
end





