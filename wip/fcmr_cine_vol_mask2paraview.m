% fcmr_4dvol_mask2paraview

%% Inputs
reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr202';
frameES  = 9;
frameED  = 21;
maskFileName = 'mask_function_LV';

%% Setup
cineVolDir   = fullfile( reconDir, 'cine_vol' );
cineVol4dDir = fullfile( reconDir, 'cine_vol_4d' );

cd(reconDir);
mkdir(cineVol4dDir);

copyfile( fullfile( cineVolDir, 'cine_vol.nii.gz'), cineVol4dDir  );


%% Load Data
cd(cineVol4dDir);

gunzip( 'cine_vol.nii.gz' );
gunzip( [maskFileName, '.nii.gz'] );

cine_vol = niftiread('cine_vol.nii');
info     = niftiread('cine_vol.nii');
nF       = size(cine_vol,4);

mask     = niftiread([maskFileName, '.nii.gz']);


%% Write .vtk files

[x,y,z] = meshgrid( 1:size(cine_vol,2), 1:size(cine_vol,1), 1:size(cine_vol,3) );

for tt = 1:nF
    
    % magnitude cine volume
    vtkwrite(['cine_vol_t' num2str(tt-1) '.vtk'], ...
             'structured_grid', y, x, z, ... 
             'scalars', 'magnitude_intensity', cine_vol(:,:,:,tt) );
         
    % mask
    vtkwrite([maskFileName '_t' num2str(tt-1) '.vtk'], ...
             'structured_grid', y, x, z, ... 
             'scalars', 'magnitude_intensity', mask(:,:,:,tt) );

end


%% Write .vtk with ES segmentation at ED timepoint

maskESatED = mask;
maskESatED(:,:,:,frameES) = 0;
maskESatED(:,:,:,frameED) = mask(:,:,:,frameES);

for tt = 1:nF
         
    % mask
    vtkwrite([maskFileName '_ESatED_t' num2str(tt-1) '.vtk'], ...
             'structured_grid', y, x, z, ... 
             'scalars', 'magnitude_intensity', maskESatED(:,:,:,tt) );

end
