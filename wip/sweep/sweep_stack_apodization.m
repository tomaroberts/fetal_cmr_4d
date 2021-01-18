function S = sweep_stack_apodization( S, varargin )
%SWEEP_STACK_APODIZATION  apodize sweep stack to remove k-t recon corrupted slices
%
%   S = SWEEP_STACK_APODIZATION( S ) 
%
%   SWEEP_STACK_APODIZATION( ..., 'name', value ) specifies optional input argument. 
%   See code for name-value pairs.
%
%   - Default apodizationlength = 0. Must be > 0 to perform apodization.

%   tar (t.roberts@kcl.ac.uk)


%% Optional Input Argument Default Values

default.reconDir        = pwd;
default.resultsDir      = pwd;
default.apodLength      = 0;
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

add_param_fn(  p, 'apodizationlength', default.apodLength, ...
        @(x) validateattributes( x, {'numeric'}, {'scalar','nonnegative'}, mfilename ) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, S, varargin{:} );

reconDir        = p.Results.recondir;
resultsDir      = p.Results.resultsdir;
apodLength      = p.Results.apodizationlength;
isVerbose       = p.Results.verbose;


%% Setup

dataDir     = fullfile( reconDir, 'data' );
maskDir     = fullfile( reconDir, 'mask' );
ktreconDir  = fullfile( reconDir, 'ktrecon' );
isVerbose   = false;

nStack      = numel( S );


%% Load Sweep Data

for iStk = 1:nStack

    % Load Params
    MAT = matfile( S(iStk).rltParamFile );
    P(iStk) = MAT.PARAM;

    % Load NIfTI
    R(iStk) = load_untouch_nii( S(iStk).rltAbFile );
    
    clear MAT
    
end
    

%% Apodize Sweep Windows
% Filter out artefacted images at start/end of sweep window
% Small side-effect: this will trim the Sweep volume by 2 * apodizationLength
% - e.g.: if Sweep volume has 1024 locations and apodizationLength = 12, then
%   the first and last 12 locations will be trimmed away
% - Cardiac phases etc. must be updated accordingly

% apodLength = 0;

for iStk = 1:nStack
    
    nX   = size( R(iStk).img, 1 );
    nY   = size( R(iStk).img, 2 );
    nZ   = size( R(iStk).img, 3 );
    nT   = size( R(iStk).img, 4 );
    
    % TODO: automate apodization based on median + SD of signal differences
    [ R(iStk).img, P(iStk).Sweep ] = sweep_window_filter( R(iStk).img, P(iStk).Sweep, apodLength );
    
    % Update Nifti - 4-D to 3-D
    % dim
    R(iStk).hdr.dime.dim(1) = 3;
    R(iStk).hdr.dime.dim(4) = size( R(iStk).img,3 );
    R(iStk).hdr.dime.dim(5) = 1;
    
    % pixdim
    R(iStk).hdr.dime.pixdim(4) = R(iStk).hdr.dime.pixdim(4) / nT;
    
    % affine - adjust z-location of first slice due to apodization
    % New zLoc(1) = Old zLoc(1) + (apodLength * z-step)
    R(iStk).hdr.hist.srow_z(4) = R(iStk).hdr.hist.srow_z(4) + ( apodLength * R(iStk).hdr.dime.pixdim(4) );
    
    % Save
    save_untouch_nii( R(iStk), [S(iStk).rltAbFile(1:end-7) '_swp3d_apod.nii.gz' ] );
    
end


%% Resample Mask to Apodized Stack Dimensions
% TODO: think nifti headers may be wrong - update to match apodized image
% volumes above

for iStk = 1:nStack

    % Load Mask
    M(iStk) = load_untouch_nii( S(iStk).maskHeartFile );
       
    % Resample
    nX   =  R(iStk).hdr.dime.dim(2);
    nY   =  R(iStk).hdr.dime.dim(3);
    nZnT =  R(iStk).hdr.dime.dim(4); 
    
    M(iStk).img = imresize3( M(iStk).img, [nX nY nZnT] );
    
    % Update Nifti Header
    M(iStk).hdr.dime.dim    = R(iStk).hdr.dime.dim;
    M(iStk).hdr.dime.pixdim = R(iStk).hdr.dime.pixdim;
    M(iStk).hdr.hist        = R(iStk).hdr.hist;
    
    % Save
    save_untouch_nii( M(iStk), [S(iStk).maskHeartFile(1:end-7) '_swp3d_apod.nii.gz' ] );   
    
end


%% Update S

for iStk = 1:nStack
    
    S(iStk).apodLength = apodLength;
    S(iStk).isApod( [1:apodLength, end-apodLength+1:end] ) = 1;
    
    S(iStk).rltAbSwpApodFile     = [S(iStk).rltAbFile(1:end-7) '_swp3d_apod.nii.gz' ];
    S(iStk).maskHeartSwpApodFile = [S(iStk).maskHeartFile(1:end-7) '_swp3d_apod.nii.gz' ];
    
end


%% Save Results to .mat File

save( fullfile( resultsDir, 'results_cardsync_sweep_intrastack.mat' ), 'S', '-v7.3' );


% sweep_stack_apodization()
end