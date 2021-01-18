

%% Register cine_vols and transform vel_vols

%%
fcmrNum = 189;

% Identify ED / ES
% - currently doing this manually in MITK
cardPhaseOrig = [16 6]; % not sure temporally identical?
cardPhaseHHH  = [9 19];


%%
fcmrNum = 191;

% Identify ED / ES
% - currently doing this manually in MITK
cardPhaseOrig = [11 24];
cardPhaseHHH  = [5 18];


%%
fcmrNum = 194;

% Identify ED / ES
% - currently doing this manually in MITK
cardPhaseOrig = [16 3];
cardPhaseHHH  = [21 8];


%%
fcmrNum = 202;

% Identify ED / ES
% - currently doing this manually in MITK
cardPhaseOrig = [20 7];
cardPhaseHHH  = [20 7];


%%
fcmrNum = 214;

% Identify ED / ES
% - currently doing this manually in MITK
cardPhaseOrig = [13 1];
cardPhaseHHH  = [13 1];


%% Load Data

reconDirOrig  = ['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum)];
reconDirHHH   = ['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum) '_hrh_fullRecon'];
cd(reconDirHHH);

outputDirName = '4dvol_reg';
outputDirPath = fullfile( reconDirHHH, outputDirName ); 
mkdir(outputDirPath);

% magnitude
cineVolDir = 'cine_vol';
cineVolFileName = 'cine_vol-RESLICE.nii.gz';

cineVolOrig = load_untouch_nii( fullfile( reconDirOrig, cineVolDir, cineVolFileName ) );
cineVolHHH  = load_untouch_nii( fullfile( reconDirHHH, cineVolDir, cineVolFileName ) );

% velocity
velVolDir = 'vel_vol';

velVolOrig0 = load_untouch_nii( fullfile( reconDirOrig, velVolDir, 'velocity-final-polyCorr-RESLICE-0.nii.gz' ) );
velVolOrig1 = load_untouch_nii( fullfile( reconDirOrig, velVolDir, 'velocity-final-polyCorr-RESLICE-1.nii.gz' ) );
velVolOrig2 = load_untouch_nii( fullfile( reconDirOrig, velVolDir, 'velocity-final-polyCorr-RESLICE-2.nii.gz' ) );
velVolHHH0 = load_untouch_nii( fullfile( reconDirHHH, velVolDir, 'velocity-final-polyCorr-RESLICE-0.nii.gz' ) );
velVolHHH1 = load_untouch_nii( fullfile( reconDirHHH, velVolDir, 'velocity-final-polyCorr-RESLICE-1.nii.gz' ) );
velVolHHH2 = load_untouch_nii( fullfile( reconDirHHH, velVolDir, 'velocity-final-polyCorr-RESLICE-2.nii.gz' ) );


%% Align Cardiac Phases

if abs(diff( cardPhaseOrig )) ~= abs(diff( cardPhaseHHH ))
    error('Must identical ED and ES separation between two volumes.');
end

% Match Cardiac Phases

nFrames = size( cineVolOrig.img, 4 );

if size( cineVolOrig.img, 4 ) ~= size( cineVolHHH.img, 4 )
    error('Number of frames not consistent between volumes.');
end

% magnitude
cineVolOrig.img = circshift( cineVolOrig.img, (nFrames - cardPhaseOrig(1)), 4 );
cineVolHHH.img  = circshift( cineVolHHH.img, (nFrames - cardPhaseHHH(1)), 4 );

% velocity
velVolOrig0.img = circshift( velVolOrig0.img, (nFrames - cardPhaseOrig(1)), 4 );
velVolOrig1.img = circshift( velVolOrig1.img, (nFrames - cardPhaseOrig(1)), 4 );
velVolOrig2.img = circshift( velVolOrig2.img, (nFrames - cardPhaseOrig(1)), 4 );
velVolHHH0.img  = circshift( velVolHHH0.img, (nFrames - cardPhaseHHH(1)), 4 );
velVolHHH1.img  = circshift( velVolHHH1.img, (nFrames - cardPhaseHHH(1)), 4 );
velVolHHH2.img  = circshift( velVolHHH2.img, (nFrames - cardPhaseHHH(1)), 4 );


%% Save Shifted Volumes
cd(outputDirPath);
save_untouch_nii( cineVolOrig, 'cine_vol_orig.nii.gz' );
save_untouch_nii( cineVolHHH, 'cine_vol_hhh.nii.gz' );

save_untouch_nii( velVolOrig0, 'vel_vol_0_orig.nii.gz' );
save_untouch_nii( velVolOrig1, 'vel_vol_1_orig.nii.gz' );
save_untouch_nii( velVolOrig2, 'vel_vol_2_orig.nii.gz' );
save_untouch_nii( velVolHHH0, 'vel_vol_0_hhh.nii.gz' );
save_untouch_nii( velVolHHH1, 'vel_vol_1_hhh.nii.gz' );
save_untouch_nii( velVolHHH2, 'vel_vol_2_hhh.nii.gz' );

mkdir roi


%% Copy across ROIs David drew before
cd(reconDirOrig);

% Folder check
if exist('roi_v2') == 7
    roiDir = 'roi_v2';
    cd(roiDir);
else
    disp('Version 2 folder does not exist ... ')
    roiDir = 'roi';
    cd(roiDir);
end

dirNames = dir('*.nii.gz');
for dd = 1:numel(dirNames)
    copyfile(dirNames(dd).name, [outputDirPath, '/roi']);
end

cd(outputDirPath);

%% Run Registration Script using Powershell

% TODO: figure out how to run from within MATLAB
% - tricky because need to call powershell to run bash script

% system('powershell -inputformat none ./register_cine_vol.bash .');
% system('powershell -inputformat none ./register_cine_vol.ps1 .');





