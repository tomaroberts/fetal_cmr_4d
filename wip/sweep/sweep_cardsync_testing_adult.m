%% sweep cardsync_testing_adult
%
% Working through 4d recon pipeline with M2D data to figure out new
% cardsync approach
%
% For now: Adult dataset, reconstructed to have 20 slices (74 frames each)
%

reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Adult\2019_03_28_Adult_Heart_m2d_w74';
reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Adult\2019_03_28_Adult_Heart';
cd(reconDir);


%% Step 0)

% k-t recon with 20 slices, each with 74 frames
% custom, upsampled CSM. Need to fix this eventually.


%% Convert sweep ktrecon output to '3.5D' (i.e.: time/space in 3rd dimension) volume

% % TODO: fold into preproc.m ?
% 
% ktreconDir  = fullfile( reconDir, 'ktrecon' );
% rltFileList = dir( fullfile( ktreconDir, 's*_rlt_ab.nii.gz' ) );
% nStack = numel(rltFileList);
% 
% for iStk = 1:nStack
% 
%     R = load_untouch_nii( fullfile( ktreconDir, rltFileList(iStk).name ) );
%     
%     nX = size( R.img, 1 );
%     nY = size( R.img, 2 );
%     n3 = size( R.img, 3 ) * size( R.img, 4 );
%     
%     SWP.img = reshape( permute(R.img,[1,2,4,3]), nX, nY, n3 );
%     
%     %TODO: save new swp *.nii equivalent to LJ python output
%     %TODO: retain *_rlt_* naming convention or rename *_swp_* ?
%     
% end
% 
% % TODO:
% % copyfile( '/ktrecon', '/data' );


%% Copy files to /data
copyfile( 'ktrecon/s10_rlt_ab.nii.gz', 'data/s10_rlt_ab.nii.gz' );
copyfile( 'ktrecon/s10_dc_ab.nii.gz',  'data/s10_dc_ab.nii.gz' );


%% Run preproc.m
cd(reconDir);
hrRange   = [40 110];
acqMethod = 'm2d'; % or: 'swp'
S = preproc( reconDir, hrRange, acqMethod );
save( fullfile( reconDir, 'data', 'results.mat' ), 'S', '-v7.3' );
disp('Preproc complete ...');


%% Run cardsync_intraslice.m
cd(reconDir);
hrRange     = [40 110];
reconDir    = pwd;
dataDir     = fullfile( reconDir, 'data' );
cardsyncDir = fullfile( reconDir, 'cardsync' );
M = matfile( fullfile( dataDir, 'results.mat' ) );
disp('Cardsync_intraslice running ...');
S = cardsync_intraslice( M.S, 'resultsDir', cardsyncDir, 'hrrange', hrRange, 'verbose', true );
disp('Cardsync_intraslice complete ...');


















