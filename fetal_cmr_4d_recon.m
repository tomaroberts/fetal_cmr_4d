%FETAL_CMR_4D_RECON  wrapper to perform fetal_cmr_4d reconstruction
%
%   Wrapper script to perform fetal whole-heart 4d volumetric reconstruction
%
%   NOTE: this will not work as a continuous script. It is designed to show
%   the reconstruction process and be run step-by-step, switching to Bash
%   when necessary for running the C++ code using SVRTK.
%
%   See the tutorial readme for more additional guidance.
%

% tar (t.roberts@kcl.ac.uk)
% jfpva (joshua.vanamerom@kcl.ac.uk)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 0 --- Perform kt SENSE reconstruction on acquired MRI data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reconDir = 'C:\fcmr_4d_recons\fcmr100';
cd(reconDir);

%% 1) Create working directories:
% In shell:
% RECONDIR=~/path/to/recon/directory
% cd $RECONDIR
% mkdir data ktrecon mask cardsync raw

% reconFrame reconstruction: copy raw data to /raw directory

%% 2) Reconstruct real-time data using reconFrame
%
% Once data is reconstructed, copy reconstructed stacks to /data
% 
% In shell:
% cp ktrecon/s*_dc_ab.nii.gz data;
% cp ktrecon/s*_rlt_ab.nii.gz data;





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1 --- 4D Magnitude CINE volume reconstruction:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3) Manually draw fetal heart masks for each stack in MITK
%       - use sXX_dc_ab.nii.gz
%       - save as sXX_mask_heart.nii.gz in /mask

%% 4) Run preproc.m
cd(reconDir);
S = preproc( reconDir );
save( fullfile( reconDir, 'data', 'results.mat' ), 'S', '-v7.3' );
disp('Preproc complete ...');

%% 5) Run recon_ref_vol.bash
% In shell:
% RECONDIR=~/path/to/recon/directory
% mkdir $RECONDIR/ref_vol
% ./recon_ref_vol.bash $RECONDIR ref_vol

%% 6) Draw fetal chest ROI using ref_vol.nii.gz
%       - save as mask_chest.nii.gz in /mask

%% 7) Run recon_dc_vol.bash
% In shell:
% RECONDIR=~/path/to/recon/directory
% ./recon_dc_vol.bash $RECONDIR dc_vol

%% 8) Run cardsync_intraslice.m
cd(reconDir);
reconDir    = pwd;
dataDir     = fullfile( reconDir, 'data' );
cardsyncDir = fullfile( reconDir, 'cardsync' );
M = matfile( fullfile( dataDir, 'results.mat' ) );
disp('Cardsync_intraslice running ...');
S = cardsync_intraslice( M.S, 'resultsDir', cardsyncDir, 'verbose', true );
disp('Cardsync_intraslice complete ...');

%% 9) Run recon_slice_cine.bash
% In shell:
% RECONDIR=~/path/to/recon/directory
% ./recon_slice_cine.bash $RECONDIR

%% 10) Run cardsync_interslice.m
%- this outputs /cardsync/results_interslice_cardsync.mat
%               /cardsync/cardphases_interslice_cardsync.txt
%               /cardsync/log_cardsync_interslice.txt

cd(reconDir);

% setup
reconDir    = pwd;
dataDir     = fullfile( reconDir, 'data' );
cardsyncDir = fullfile( reconDir, 'cardsync' );
cineDir     = fullfile( reconDir, 'slice_cine_vol' );    
M = matfile( fullfile( cardsyncDir, 'results_cardsync_intraslice.mat' ) );

% target slice
tgtLoc = NaN;
tgtLocFile = fullfile( dataDir, 'tgt_slice_no.txt' );
if exist( tgtLocFile , 'file' )
  fid = fopen( tgtLocFile, 'r' );
  tgtLoc = fscanf( fid, '%f' );
  fclose( fid );
end

% excluded slices
excludeSlice = [];
excludeSliceFile = fullfile( dataDir, 'force_exclude_slice.txt' );
if exist( excludeSliceFile , 'file' )
  fid = fopen( excludeSliceFile, 'r' );
  excludeSlice = fscanf( fid, '%f' ) + 1;  % NOTE: slice locations in input file are zero-indexed
  fclose( fid );
end

% slice-slice cardiac synchronisation
S = cardsync_interslice( M.S, 'recondir', cineDir, 'resultsdir', cardsyncDir, 'tgtloc', tgtLoc, 'excludeloc', excludeSlice );
disp('Cardsync_interslice complete ...');

%% 11) Run recon_cine_vol.bash
% In shell:
% RECONDIR=~/path/to/recon/directory
% ./recon_cine_vol.bash $RECONDIR cine_vol

%% 12) Summarise in MATLAB
% (non-essential)
S = summarise_recon( [reconDir '/cine_vol'], [reconDir '/cardsync'], 'verbose', true );
I = plot_info( [reconDir '/cine_vol/info.tsv'] );





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2 --- 4D Velocity CINE volume reconstruction:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 13) Manually draw uterus masks for each stack in MITK
%       - use sXX_dc_ab.nii.gz
%       - save as sXX_mask_uterus.nii.gz in /mask

%% 14) Run fcmr_4dflow_preprocessing.m
cd(reconDir);
fcmr_4dflow_preprocessing( reconDir );
disp('fcmr_4dflow_preprocessing complete ...');

%% 15) Run fcmr_4dflow_get_first_moments.m

% - NOTE: this script creates gradient first moment .txt files using 
% data created by reconFrame. Alternatively, the two .txt files required 
% for velocity volume reconstruction can be calculated manually if you know
% the parameters of your flow encoding gradients

cd(reconDir);
fcmr_4dflow_get_first_moments( reconDir );
disp('fcmr_4dflow_get_first_moments complete ...');

%% 16) Run recon_vel_vol.bash
% In shell:
% RECONDIR=~/path/to/recon/directory
% ./recon_vel_vol.bash $RECONDIR vel_vol

%% 17) Run fcmr_4dflow_postprocessing.m

% - TODO: update script with MATLAB niftiread/niftiwrite to avoid having to
% use reslice_nii.m

cd(reconDir);
bloodpoolMask = 'mask_blood_pool';
velMasks      = bloodpoolMask;
fcmr_4dflow_postprocessing( reconDir, 'useVelDriftCorr', true, 'fileExt', 'polyCorr', 'bloodpoolMask', bloodpoolMask, 'velMasks', velMasks );
disp('fcmr_4dflow_postprocessing complete ...');




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3 --- Whole-heart 4D visualisation in Paraview:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 18) Run fcmr_4dflow_make_vector_vol.py
%
% In python:
% - edit fcmr_4dflow_make_vector_vol.py with correct path information
% Run: fcmr_4dflow_make_vector_vol.py

%% 19) Load state file in Paraview.
% /vel_vol_4d/paraview_polyCorr contains .pvsm file
% In Paraview:
% - File > Load State > fcmr***_paraview.pvsm



