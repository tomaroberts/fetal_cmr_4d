%% sweep cardsync_testing_adult
%
% Working through 4d recon pipeline with M2D data to figure out new
% cardsync approach
%
% For now: Adult dataset, 5 stacks, reconstructed to have 20 slices (74 frames each)
%

seriesNos = [11 14 16 18 20];

reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Adult\2019_03_28_Adult_Heart_swp_w74';
cd(reconDir);

mkdir data
mkdir cardsync
mkdir mask


%% Step 0)

% k-t recon with 20 slices, each with 74 frames
% custom, upsampled CSM. Need to fix this eventually.


%% Convert sweep ktrecon output to '3.5D' (i.e.: time/space in 3rd dimension) volume

% TODO: fold into preproc.m ?

cd([reconDir '/ktrecon']);
for seriesNo = seriesNos
    seriesNoStr = num2str(seriesNo);    
    swp_unfold( ['s' seriesNoStr '_rlt_ab.nii.gz'], ['s' seriesNoStr '_rlt_ab_swp3d.nii.gz']);
end


%% Copy files to /data
cd(reconDir);
for seriesNo = seriesNos
    seriesNoStr = num2str(seriesNo); 
    copyfile( ['ktrecon/s' seriesNoStr '_rlt_ab.nii.gz'], ['data/s' seriesNoStr '_rlt_ab.nii.gz'] );
    copyfile( ['ktrecon/s' seriesNoStr '_dc_ab.nii.gz'],  ['data/s' seriesNoStr '_dc_ab.nii.gz'] );
end


%% Run preproc.m
cd(reconDir);
hrRange   = [65 85];
acqMethod = 'm2d'; % or: 'swp'
S = preproc( reconDir, hrRange, acqMethod );
close all;
save( fullfile( reconDir, 'data', 'results.mat' ), 'S', '-v7.3' );
disp('Preproc complete ...');


%% Run cardsync_intraslice.m
cd(reconDir);
hrRange     = [65 85];
reconDir    = pwd;
dataDir     = fullfile( reconDir, 'data' );
cardsyncDir = fullfile( reconDir, 'cardsync' );
M = matfile( fullfile( dataDir, 'results.mat' ) );
disp('Cardsync_intraslice running ...');
S = cardsync_intraslice( M.S, 'resultsDir', cardsyncDir, 'hrrange', hrRange, 'verbose', true );
disp('Cardsync_intraslice complete ...');


















