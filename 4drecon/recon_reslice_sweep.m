function S = recon_reslice_sweep( S, varargin )
%RECON_RESLICE_SWEEP  reslice sweep volume for SVRTK
%
%   S = RECON_RESLICE_SWEEP( S ) 
%
%   RECON_RESLICE_SWEEP( ..., 'name', value ) specifies optional input argument. 
%   See code for name-value pairs.
%

%   tar (t.roberts@kcl.ac.uk)


%% Optional Input Argument Default Values

default.reconDir        = pwd;
default.resultsDir      = pwd;
default.nDyn            = 64;
default.isVerbose       = true;


%% Parse Input

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end


addRequired(   p, 'S', ... 
    @(x) validateattributes( x, {'struct'}, {'vector'}, mfilename) );

add_param_fn(   p, 'recondir', default.reconDir, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(   p, 'resultsdir', default.resultsDir, ...
    @(x) validateattributes( x, {'char'}, {'vector'}, mfilename) );

add_param_fn(  p, 'ndyn', default.nDyn, ...
        @(x) validateattributes( x, {'numeric'}, {'scalar','nonnegative'}, mfilename ) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, S, varargin{:} );

reconDir        = p.Results.recondir;
resultsDir      = p.Results.resultsdir;
nDyn            = p.Results.ndyn;
isVerbose       = p.Results.verbose;


%% Setup

dataDir     = fullfile( reconDir, 'data' );
maskDir     = fullfile( reconDir, 'mask' );
ktreconDir  = fullfile( reconDir, 'ktrecon' );
isVerbose   = false;

nStack      = numel( S );


%% Load Sweep Stacks and Masks

for iStk = 1:nStack

    % Load Params
    MAT = matfile( S(iStk).rltParamFile );
    P(iStk) = MAT.PARAM;

    % Load Sweep Data
    R(iStk) = load_untouch_nii( S(iStk).rltAbSwpApodFile );
    
    % Load Masks
    M(iStk) = load_untouch_nii( S(iStk).maskHeartSwpApodFile );
    M(iStk).img = single( M(iStk).img );

    clear N MAT
    
end


%% Reslice Sweep Volumes into 4-D

for iStk = 1:nStack 
    
    % Reslice Configuration
    numSwpLoc = size( R(iStk).img, 3 );
    nSlices   = floor( numSwpLoc / nDyn ); % floor - nSlices must be integer
    
    % Trim Sweep Stack to allow discretized slices
    trimLength(iStk) = ( numSwpLoc - ( nSlices * nDyn ) ) / 2;
    S(iStk).isTrimSwpLoc = zeros( size(S(iStk).thetaFrameSwpLoc) );
    S(iStk).isTrimSwpLoc( [1:trimLength(iStk), end-trimLength(iStk)+1:end] ) = 1;
    
    R(iStk).img( :,:,S(iStk).isTrimSwpLoc==1 ) = [];
    M(iStk).img( :,:,S(iStk).isTrimSwpLoc==1 ) = [];
    
    % Redefine after Trim
    numSwpLoc = size( R(iStk).img, 3 );
    binWidthSlices = numSwpLoc / nSlices;
    
    % Init Slice Bins
    edgesSlices = 1:binWidthSlices:numSwpLoc+binWidthSlices;
    
    % Reslice Sweep Volume   
    for iS = 1:nSlices

        currentBinRange = edgesSlices(iS):edgesSlices(iS+1)-1;
        
        R_resliced(iStk).img( :,:,iS,: ) = R(iStk).img( :,:,currentBinRange );
        M_resliced(iStk).img( :,:,iS,: ) = M(iStk).img( :,:,currentBinRange );

    end
    
end


%% Save NIfTI

for iStk = 1:nStack
    
    % Sweep Stack
    R(iStk).img = R_resliced(iStk).img;
    
    % Update Header
    nZ = size( R_resliced(iStk).img, 3 );
    nT = size( R_resliced(iStk).img, 4 );
    
    % dim
    R(iStk).hdr.dime.dim(1)     = 4;
    R(iStk).hdr.dime.dim([4,5]) = [nZ, nT];
       
    % affine - adjust z-location of first slice due any trimming
    % New zLoc(1) = Old zLoc(1) + (trimLength * z-step)
    R(iStk).hdr.hist.srow_z(4) = R(iStk).hdr.hist.srow_z(4) + ( trimLength(iStk) * R(iStk).hdr.dime.pixdim(4) );

    % pixdim
    warning('REMINDER: R.hdr.pixdim(4) hard-coded to 4. This is probably very wrong.');
    R(iStk).hdr.dime.pixdim(4) = 4; % Standard slice thk with 2mm overlap
%     R(iStk).hdr.dime.pixdim(4) = R(iStk).hdr.dime.pixdim(4) * nT;   
    
    % Save Resliced Sweep
    S(iStk).rltAbSwpReslicedFile = fullfile( dataDir, sprintf( '%s_rlt_ab_swp_resliced.nii.gz', S(iStk).desc ) );
    save_untouch_nii( R(iStk), S(iStk).rltAbSwpReslicedFile );
    
    
    
    % Mask
    % NB: as resampling from full, 3-D mask, we have every dynamic masked
    % in the resliced 4-D nifti. This is different to M2D maskHeartFile,
    % which only has first frames masked.
    
    % 4-D Mask
%     M(iStk).img = M_resliced(iStk).img;
%     
%     % Update Header
%     M(iStk).hdr.dime.dim    = R(iStk).hdr.dime.dim;
%     M(iStk).hdr.dime.pixdim = R(iStk).hdr.dime.pixdim;
%     M(iStk).hdr.hist        = R(iStk).hdr.hist;

    % 3-D Mask
    M(iStk).img = M_resliced(iStk).img( :,:,:,round(nDyn/2) );
    
    % Update Header
    M(iStk).hdr.dime.dim([1, 4]) = [3, nZ];
    warning('REMINDER: M.hdr.pixdim(4) hard-coded to 4. This is probably very wrong.');
    M(iStk).hdr.dime.pixdim(4)   = 4; % Standard slice thk with 2mm overlap
    M(iStk).hdr.dime.pixdim(5)   = 1;
    M(iStk).hdr.hist             = R(iStk).hdr.hist;
    
    % Save
    S(iStk).maskHeartSwpReslicedFile = fullfile( maskDir, sprintf( '%s_mask_heart_swp_resliced.nii.gz', S(iStk).desc ) );
    save_untouch_nii( M(iStk), S(iStk).maskHeartSwpReslicedFile );
    
end


%% Update thetaFrame Timings
for iStk = 1:nStack
    
    S(iStk).thetaFrameSwpLocResliced = S(iStk).thetaFrameSwpLoc;
    S(iStk).thetaFrameSwpLocResliced( S(iStk).isTrimSwpLoc==1 ) = [];
    
end


%% Resample tRR
% Think this seems reasonable...

for iStk = 1:nStack

    S(iStk).hrMeanSwpResliced = S(iStk).hrMeanSwpLoc;
    
    % Remove Apodized Frames
    S(iStk).hrMeanSwpResliced( S(iStk).isApod==1 ) = [];
    
    % Remove Trimmed Frames
    S(iStk).hrMeanSwpResliced( S(iStk).isTrimSwpLoc==1 ) = [];
    
    % Calculate tRR per Slice
    S(iStk).tRRSwpResliced = 60 ./ mean( reshape( S(iStk).hrMeanSwpResliced, nDyn, [] ), 1);    
    
end


%% Save Results to Text Files

fid = fopen( fullfile( resultsDir, 'mean_rrinterval.txt' ), 'w' );
fprintf( fid, '%.6f ', mean( cat( 2, S.tRRSwpResliced ) ) );
fclose( fid );
fid = fopen( fullfile( resultsDir, 'rrintervals.txt' ), 'w' );
fprintf( fid, '%.6f ', cat( 2, S.tRRSwpResliced ) );
fclose( fid );
fid = fopen( fullfile( resultsDir, 'cardphases_swp_resliced_cardsync.txt' ), 'w' );
fprintf( fid, '%.6f ', cat( 1, S.thetaFrameSwpLocResliced ) );
fclose( fid );


%% Save Results to .mat File

save( fullfile( resultsDir, 'results_cardsync_sweep_resliced.mat' ), 'S', '-v7.3' );



% recon_reslice_sweep(...)
end
