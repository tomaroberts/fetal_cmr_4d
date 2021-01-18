%FETAL_CMR_4D_RECON_ADAPTIVE  wrapper to perform adaptive fetal_cmr_4d reconstruction
%
%   Wrapper script to perform fetal whole-heart 4d volumetric
%   reconstruction with adaptive regularization
%
%   NOTE: this will not work as a continuous script. It is designed to show
%   the reconstruction process and be run step-by-step, switching to Bash
%   when necessary for running the C++ code using SVRTK.
%
%   See the tutorial readme for more additional guidance.
%

% tar (t.roberts@kcl.ac.uk)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1 --- Perform uniform regularization kt SENSE reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fcmrNum = 189;

reconDir    = ['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum)];
reconHrhDir = [ reconDir, '_hrh' ];

cd(reconDir);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2 --- Perform hierarchical kt SENSE reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ArrangeKT on beastie01
% SolverKT on gpubeastie03
% recon_exam_hrh on beastie01


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3 --- Perform harmonic regularization kt SENSE reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% recon_exam using xtTrnCus / harmonic regularization on beastie01


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4 --- Perform 4D CINE volume reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3) Copy s*_adaptive.nii.gz files to /data

cd(reconHrhDir);

mkdir data
mkdir cine_vol_harmonic

dcFiles  = dir('ktrecon_harmonic/s*_dc_ab.nii.gz');
rltFiles = dir('ktrecon_harmonic/s*_rlt_ab.nii.gz');

for ii = 1:numel(dcFiles)
    copyfile( ['ktrecon_harmonic/' dcFiles(ii).name], reconHrhDir );
    movefile( dcFiles(ii).name , ['data/' dcFiles(ii).name(1:end-7) '_hrm.nii.gz' ] );

    copyfile( ['ktrecon_harmonic/' rltFiles(ii).name], reconHrhDir );
    movefile( rltFiles(ii).name , ['data/' rltFiles(ii).name(1:end-7) '_hrm.nii.gz' ] );
end


%% 4) Inspect s*_ab_rlt_hrm.nii.gz stacks for alias artefacts
%TODO: remove alias artefacts entirely

% Create reconHrhDir/data/force_exclude_slice.txt


%% 5) Get excluded frames from original cine_vol reconstruction

reconHrhDir = [ reconDir, '_hrh' ];
get_cine_vol_force_exclude( reconDir );
movefile( 'force_exclude_cine_vol.txt' , [reconHrhDir '/data' ] );
cd( [reconHrhDir '/data' ] );


%% 5) Run recon_cine_vol_harmonic.bash

% Copy files to gpubeastie03:
% fcmrXXX:
%   - /cardsync
%   - /cine_vol
%   - /data
%   - /dc_vol
%
% fcmrXXX_hrh:
%   - data
%       - important: should contain force_exclude_cine_vol.txt (and
%                    force_exclude_slice.txt)
%   - cine_vol_harmonic (this will be created)

% In shell:
% RECONDIR=~/path/to/recon/directory --- ie: ~fcmrXXX_hrh
% RECONTGTDIR=~/path/to/first_recon/directory --- ie: ~/fcmrXXX
% ./recon_cine_vol_harmonic.bash $RECONDIR $RECONTGTDIR cine_vol_harmonic
%
% - this takes the transformations from the cine_vol reconstruction and
% applies them to the harmonic kt reconstructed stacks



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3 --- 4D Velocity CINE volume reconstruction:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 13) Manually draw uterus masks for each stack in MITK
%       - use sXX_dc_ab.nii.gz
%       - save as sXX_mask_uterus.nii.gz in /mask

%% 14) Run fcmr_4dflow_preprocessing.m
cd(reconHrhDir);
ktreconDir = 'ktrecon_harmonic';
cinevolDir = 'cine_vol_harmonic';
fcmr_4dflow_preprocessing( reconHrhDir, ktreconDir, cinevolDir );
disp('fcmr_4dflow_preprocessing complete ...');

%% 15) Run fcmr_4dflow_get_first_moments.m

%%% Alternative to this cell: get .txt files from existing recon

% - NOTE: this script creates gradient first moment .txt files using 
% data created by reconFrame. Alternatively, the two .txt files required 
% for velocity volume reconstruction can be calculated manually if you know
% the parameters of your flow encoding gradients

cd(reconHrhDir);
ktreconDir = 'ktrecon_harmonic';
fcmr_4dflow_get_first_moments( reconHrhDir, 'ktreconDir', ktreconDir );
disp('fcmr_4dflow_get_first_moments complete ...');

%% 16) Run recon_vel_vol.bash

% Copy reconHrhDir/data to gpubeastie03

% In shell:
% RECONDIR=~/path/to/recon/directory
% ./recon_vel_vol.bash $RECONDIR vel_vol











