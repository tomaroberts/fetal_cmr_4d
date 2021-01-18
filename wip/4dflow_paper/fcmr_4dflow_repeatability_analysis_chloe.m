%% fcmr_4dflow_repeatability_analysis --- Chloe
% - August 2020 - analysis of Chloe's ROIs
% - nb: blindingDir from paper analysis by Milou/David
%


%% Directories
studyDir = strcat('I:\King', '''', 's College London\Perinatal CMR - Documents\General\Data\fcmr_4dflow_paper');
blindingDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\#4dflow_paper_analysis_info\flow_analysis_repeatability';


%% Unblind
cd(blindingDir);
load('fcmrIDs_for_unblinding.mat');

paperIDs     = fcmrIDs([1,2,3,4,5,6,7],1); 
numCases     = numel(paperIDs);
fetalWeights = [2.11, 1.82, 0.88, 1.48, 1.82, 1.99, 1.01]; % kg, taken from 4D_Flow_Paper_Fetal_volumes_from_Milou.xls

roiDirPath = @( ii, roiDirName ) fullfile( studyDir, ['fcmr' num2str(paperIDs(ii))], roiDirName );

cd(studyDir);

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
%% Chloe
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

load( fullfile( blindingDir, 'chloe_flow_analysis.mat' ), 'CHLOE' );

%% Mean Flows

%%% IMPORTANT: order here is different to Excel sheet!
%%% Copied these values into:
%%% fcmr_4DFlow_ROI_reliability_analysis_FINAL.xlsx
vesselNames = CHLOE.fcmr189.T1.Properties.VariableNames(2:6);

% Milou
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








