%% FCMR kt-ml --- pixel_array size reduction

% - data too large for Chen/Gavin networks
% - reduce datasets in RO direction
%
% Tom Roberts (t.roberts@kcl.ac.uk)

fcmrNum = 191;
reconDir = ['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum)];
cd(reconDir);

outputDirPath = fullfile( 'I:\OneDrive - King''s College London\kt_ml', ['fcmr' num2str(fcmrNum)] );
% outputDirPath = fullfile( reconDir, 'kt_ml' );

mkdir( outputDirPath );


%% Init

isCompressed = false;

dataDir      = fullfile( reconDir, 'data' );
maskDir      = fullfile( reconDir, 'mask' );
ktreconDir   = fullfile( reconDir, 'ktrecon' );
kspaceDir    = fullfile( outputDirPath );
ktFactor     = 8;

fprintf( 'Performing FCMR ML Preprocessing ...\n' );


%% Get Mask Centroid/BoundingBox Info for Each Slice

% Identify Dynamic MR Image Series
rltFileList       = dir( fullfile( dataDir, '*_rlt_ab.nii.gz' ) );

% Get Number of Stacks
nStack = numel(rltFileList);
% nStack = 1;

% Initialise Stack Struct
S = struct([]);

% Read Data Into Stack Struct
for iStk = 1:nStack
  
    % Identify Files
    S(iStk).desc          = strrep( rltFileList(iStk).name, '_rlt_ab.nii.gz', '' );
    S(iStk).rltAbFile     = fullfile( rltFileList(iStk).folder,  rltFileList(iStk).name  );
    S(iStk).rltReFile     = fullfile( ktreconDir, strrep( rltFileList(iStk).name, 'ab', 're' ) );
    S(iStk).rltImFile     = fullfile( ktreconDir, strrep( rltFileList(iStk).name, 'ab', 'im' ) );
    S(iStk).rltMatFile    = fullfile( ktreconDir, sprintf( '%s_rlt_recon.mat', S(iStk).desc ) );
    S(iStk).rltParamFile  = fullfile( ktreconDir, sprintf( '%s_rlt_parameters.mat', S(iStk).desc ) );
    S(iStk).dcAbFile      = fullfile( dataDir, sprintf( '%s_dc_ab.nii.gz', S(iStk).desc ) );
    S(iStk).slwAbFile     = fullfile( ktreconDir, sprintf( '%s_slw_ab.nii.gz', S(iStk).desc ) );
    S(iStk).trnAbFile     = fullfile( ktreconDir, sprintf( '%s_trn_ab.nii.gz', S(iStk).desc ) );
    S(iStk).maskHeartFile = fullfile( maskDir, sprintf( '%s_mask_heart.nii.gz', S(iStk).desc ) );
    S(iStk).kspaceMatFile = fullfile( kspaceDir, sprintf( '%s_kspace.mat', S(iStk).desc ) );
    S(iStk).csmMatFile    = fullfile( ktreconDir, sprintf( '%s_csm.mat', S(iStk).desc ) );
    
    % Useful Conversion Functions
    PARAMS.dimX = 1; PARAMS.dimY = 2; PARAMS.dimT = 3; PARAMS.dimC = 4; PARAMS.dimZ = 5;
    kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, PARAMS.dimX ), PARAMS.dimY ) );
    xt2kt = @( kt ) ifftshift( ifftshift( fft2( kt ), PARAMS.dimX ), PARAMS.dimY );
    xt2xf = @( xt ) fftshift( fft( xt, [], PARAMS.dimT ), PARAMS.dimT );
    xf2xt = @( xf ) ifft( ifftshift( xf, PARAMS.dimT ), [], PARAMS.dimT );
    phase_correct     = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [+1;-1], size(k,1)/2, 1 ) ) ) ); 
    inv_phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [-1;+1], size(k,1)/2, 1 ) ) ) );
    
    % Load Parameters
    M = matfile( S(iStk).rltParamFile );
    P = M.PARAM;

    % Extract Parameters
    S(iStk).nLoc             = P.Timing.numLoc;
    S(iStk).sliceThickness   = P.Scan.RecVoxelSize(3);    
    
    fprintf( 'Getting centroid/bounding box of stack %s ...\n', S(iStk).desc );
    
    % Load CSM to get pixel array dimensions
    ktmlMatFilePath = fullfile( outputDirPath, strcat( S(iStk).desc, '_kt_ml_data.mat' ) );
    load( ktmlMatFilePath, 'csm' );
    nX = size( csm, 1 );
    nY = size( csm, 2 );
    
    % Load Heart Mask
    MH(iStk) = load_untouch_nii( S(iStk).maskHeartFile );

    % Pad mask to match kt_ml data arrays
    padSize = round( ( [nX,nY] - size(MH(iStk).img(:,:,1)) ) / 2 );
    MH(iStk).img  = single( padarray( MH(iStk).img, padSize, 0 ) );
    
%     % Sanity check mask matches xtRcn
% 	  load( ktmlMatFilePath, 'ktRcn' );
%     xtRcn = kt2xt( ktRcn );    
%     imtar( abs(xtRcn(:,:,1,1,6)) .* abs(MH.img(:,:,6)) );

    % Find centroid of mask
    for iSlice = 1:size( MH(iStk).img,3 )
        stats = regionprops( MH(iStk).img(:,:,iSlice) );
        H(iStk).Centroid(iSlice,:)    = round( stats.Centroid );
        H(iStk).BoundingBox(iSlice,:) = stats.BoundingBox([3,4]);
    end
    
    H(iStk).meanRoCentroid = round( mean( H(iStk).Centroid(:,2) ) );
    H(iStk).meanPeCentroid = round( mean( H(iStk).Centroid(:,1) ) );
    H(iStk).maxRoWidth = max(H(iStk).BoundingBox(:,2));
    H(iStk).maxPeWidth = max(H(iStk).BoundingBox(:,1));
    
end
  

%% Reduce Data Size

% size to resize to
roWidth = 50;
peWidth = 8*6;

meanRoCentroids = [H.meanRoCentroid];
meanPeCentroids = [H.meanPeCentroid];

maxRoWidth_All = max([H.maxRoWidth]);
maxPeWidth_All = max([H.maxPeWidth]);

if maxRoWidth_All > roWidth
    error('Heart is really large! So far restricting RO pixel dimension to 50 voxels. Might need to make bigger.');
end

% Readout Dimension Range for each Stack
rRo(1,:) = meanRoCentroids - roWidth/2;
rRo(2,:) = meanRoCentroids + roWidth/2 - 1;

rPe(1,:) = meanPeCentroids - peWidth/2;
rPe(2,:) = meanPeCentroids + peWidth/2 - 1;

% Resize data
for iStk = 1:nStack
    
    ktmlMatFilePath = fullfile( outputDirPath, strcat( S(iStk).desc, '_kt_ml_data.mat' ) );
    load( ktmlMatFilePath, 'ktRcn' );
    
    ktRcn = ktRcn( :,rPe(1):rPe(2),:,:,: );
    
end


    
    
% %     % Load NIfTI
% %     R = load_untouch_nii( S(iStk).rltAbFile );
% %     S(iStk).niiHdr = R.hdr;
% % 
% %     % Separate slices
% %     for iLoc = 1:S(iStk).nLoc
% % 
% %         % Dynamic Image Series
% %         S(iStk).frameDuration   = P.Timing.frameDuration;
% %         S(iStk).tFrame{iLoc}    = P.Timing.sliceTime(iLoc) + S(iStk).frameDuration * (0:(P.Encoding.NrDyn(1)-1));
% %         
% %     end
%     
%     
%     %%% Data To Save for ML
%     
%     % Reference - k-t SENSE kt8 Recon
% %     load( S(iStk).rltMatFile );
% %     xtRcn = squeeze( xtRcn );
% %     
% %     % Zero-Pad to Match ktAcq
% %     nX = size( xtRcn, 1 ) * 2; nY = size( xtRcn, 2 );
% %     padSize = round(  ( [nX,nY] - size(xtRcn(:,:,1)) ) / 2 );
% %     xtRcn   = padarray( xtRcn, padSize, 0 );
% 
%     ktRcnFileList = dir( fullfile( ktreconDir, sprintf( '%s_slice*.mat', S(iStk).desc ) ) );
%     for iSlice = 1:S(iStk).nLoc
%         load( fullfile( ktreconDir, ktRcnFileList(iSlice).name ), 'ktRcn' );
%         ktRcn = phase_correct( ktRcn );
%         ktRcnAll(:,:,:,iSlice) = ktRcn;
%     end    
%     ktRcn = permute( ktRcnAll, [1,2,3,5,4] ); % make 5-D for consistency
%     clear ktRcnAll;
%     
%     % Sensitivity Maps
%     load( S(iStk).csmMatFile , 'csm' );
%     
%     % K-Space Data
%     load( S(iStk).kspaceMatFile , 'ktAcq', 'ktTrn' );
%     ktAcq = phase_correct( ktAcq );
%     ktTrn = phase_correct( ktTrn );
%     
%     if isCompressed
%         [ktAcq, ktSmp] = compressKT( ktAcq );
%         [ktTrn, ~]     = compressKT( ktTrn );
%     else
%         ktSmp = single( sum( sum( ktAcq, PARAMS.dimC ), PARAMS.dimX ) ~= 0 );
%     end
%         
%     % Important Parameters
%     PARAMS.nFE      = size( ktAcq,1 );
%     PARAMS.nPE      = size( ktAcq,2 );
%     PARAMS.nFrames  = size( ktAcq,3 );
%     PARAMS.nCoils   = size( ktAcq,4 );
%     PARAMS.nSlices  = size( ktAcq,5 );
%     
%     % Save Matfile
%     fprintf( 'Saving stack %s ...\n', S(iStk).desc );
%     
%     ktmlMatFilePath = fullfile( outputDirPath, strcat( S(iStk).desc, '_kt_ml_data.mat' ) );
%     save( ktmlMatFilePath, 'ktRcn', 'csm', 'ktAcq', 'ktSmp', 'ktTrn', 'isCompressed', ...
%           'PARAMS', 'kt2xt', 'xt2kt', 'xt2xf', 'xf2xt', 'phase_correct', 'inv_phase_correct', ...
%           '-v7.3' );
%     
%     
%     % Clear Variables
%     clear R ktRcn csm ktAcq ktSmp ktTrn PARAMS
%     
%     fprintf( 'Saved stack %s ...\n', S(iStk).desc );
%     
% end

% fprintf( 'All stacks for FCMR %i saved. \n', fcmrNum );






