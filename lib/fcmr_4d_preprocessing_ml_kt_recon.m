%% FCMR Preprocessing for ML-based k-t reconstruction

% Required:
% - K-space data:     [nFE, nPE, nFrames, nCoils, nSlices] --- maybe need to regenerate 
% - Sensitivity maps: [nFE, nPE, 1, nCoils, nSlices] --- ktrecon/*csm.mat
% - Reference:        [nFE, nPE, nFrames, 1, nSlices] --- *rlt_recon.mat
% - Number of frequency encoding lines: nFE
% - Number of phase encoding lines:     nPE
% - Number of cardiac frames:           nFrames
% - Number of coils:                    nCoils
% - Number of slices:                   nSlices
% - Sliding window x-f priors: ?
%
%
% Tom Roberts (t.roberts@kcl.ac.uk)

fcmrNum = 214;
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


%% Identify Data and Get Parameters

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
    
    fprintf( 'Preprocessing stack %s ...\n', S(iStk).desc );
    
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
    
%     % Load NIfTI
%     R = load_untouch_nii( S(iStk).rltAbFile );
%     S(iStk).niiHdr = R.hdr;
% 
%     % Separate slices
%     for iLoc = 1:S(iStk).nLoc
% 
%         % Dynamic Image Series
%         S(iStk).frameDuration   = P.Timing.frameDuration;
%         S(iStk).tFrame{iLoc}    = P.Timing.sliceTime(iLoc) + S(iStk).frameDuration * (0:(P.Encoding.NrDyn(1)-1));
%         
%     end
    
    
    %%% Data To Save for ML
    
    % Reference - k-t SENSE kt8 Recon
%     load( S(iStk).rltMatFile );
%     xtRcn = squeeze( xtRcn );
%     
%     % Zero-Pad to Match ktAcq
%     nX = size( xtRcn, 1 ) * 2; nY = size( xtRcn, 2 );
%     padSize = round(  ( [nX,nY] - size(xtRcn(:,:,1)) ) / 2 );
%     xtRcn   = padarray( xtRcn, padSize, 0 );

    ktRcnFileList = dir( fullfile( ktreconDir, sprintf( '%s_slice*.mat', S(iStk).desc ) ) );
    for iSlice = 1:S(iStk).nLoc
        load( fullfile( ktreconDir, ktRcnFileList(iSlice).name ), 'ktRcn' );
        ktRcn = phase_correct( ktRcn );
        ktRcnAll(:,:,:,iSlice) = ktRcn;
    end    
    ktRcn = permute( ktRcnAll, [1,2,3,5,4] ); % make 5-D for consistency
    clear ktRcnAll;
    
    % Sensitivity Maps
    load( S(iStk).csmMatFile , 'csm' );
    
    % K-Space Data
    load( S(iStk).kspaceMatFile , 'ktAcq', 'ktTrn' );
    ktAcq = phase_correct( ktAcq );
    ktTrn = phase_correct( ktTrn );
    
    if isCompressed
        [ktAcq, ktSmp] = compressKT( ktAcq );
        [ktTrn, ~]     = compressKT( ktTrn );
    else
        ktSmp = single( sum( sum( ktAcq, PARAMS.dimC ), PARAMS.dimX ) ~= 0 );
    end
        
    % Important Parameters
    PARAMS.nFE      = size( ktAcq,1 );
    PARAMS.nPE      = size( ktAcq,2 );
    PARAMS.nFrames  = size( ktAcq,3 );
    PARAMS.nCoils   = size( ktAcq,4 );
    PARAMS.nSlices  = size( ktAcq,5 );
    
    % Save Matfile
    fprintf( 'Saving stack %s ...\n', S(iStk).desc );
    
    ktmlMatFilePath = fullfile( outputDirPath, strcat( S(iStk).desc, '_kt_ml_data.mat' ) );
    save( ktmlMatFilePath, 'ktRcn', 'csm', 'ktAcq', 'ktSmp', 'ktTrn', 'isCompressed', ...
          'PARAMS', 'kt2xt', 'xt2kt', 'xt2xf', 'xf2xt', 'phase_correct', 'inv_phase_correct', ...
          '-v7.3' );
    
    
    % Clear Variables
    clear R ktRcn csm ktAcq ktSmp ktTrn PARAMS
    
    fprintf( 'Saved stack %s ...\n', S(iStk).desc );
    
end

fprintf( 'All stacks for FCMR %i saved. \n', fcmrNum );






