%% 2019_07_21 --- QFlow --- FlowPhantom_9 --- spherical flask
%
% QFlow from .nii
% 1.25mm isotropic voxels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Step 1)
%% separate magnitude and magnitude2 images in ISDPACS nifti files

% reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow';
% cd(reconDir);
% cd data
% 
% nii = load_untouch_nii('*FH*it0.nii.gz');
% nii.img = nii.img(:,:,1:2:end);
% nii.hdr.dime.dim(4) = size(nii.img,3);
% save_untouch_nii(nii,'fh_mag.nii.gz');
% 
% nii = load_untouch_nii('*RL*it0.nii.gz');
% nii.img = nii.img(:,:,1:2:end);
% nii.hdr.dime.dim(4) = size(nii.img,3);
% save_untouch_nii(nii,'rl_mag.nii.gz');
% 
% nii = load_untouch_nii('*AP*it0.nii.gz');
% nii.img = nii.img(:,:,1:2:end);
% nii.hdr.dime.dim(4) = size(nii.img,3);
% save_untouch_nii(nii,'ap_mag.nii.gz');
% 
% % *it3.nii.gz files renamed as fh_vel.nii.gz, etc.


%% Step 2)
%% Load stacks and scale
reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow';
cd(reconDir);
cd data_ffd

magFiles = dir('*_mag*.nii.gz');
magFiles = {magFiles([1,2,3]).name}; %re-order: ap/fh/rl

velFiles = dir('*_vel*.nii.gz');
velFiles = {velFiles([1,2,3]).name};

for ii = 1:numel(magFiles)
    nii = load_untouch_nii(magFiles{ii});
    MAG(:,:,:,ii) = double(nii.img);
    
    nii = load_untouch_nii(velFiles{ii});
    VEL(:,:,:,ii) = double(nii.img);
    
    % taken from nii header:
    slope = 0.0147;
    inter = -30;
    VEL(:,:,:,ii) = ( VEL(:,:,:,ii) .* slope ) + inter;
    
    nii.img = VEL(:,:,:,ii);
    nii.hdr.dime.datatype = 16;
    save_untouch_nii(nii,[velFiles{ii}(1:end-7) '_scaled.nii.gz']);
    
end
clear nii


implay_RR([MAG(:,:,:,1), ...
           MAG(:,:,:,2), ...
           MAG(:,:,:,3)] ...
           ,'gray');

implay_RR([VEL(:,:,:,1), ...
           VEL(:,:,:,2), ...
           VEL(:,:,:,3)] ...
           ,'jet',[-5,5]);

       
       
