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

add_param_fn(  p, 'ndyn', default.nDyn, ...
        @(x) validateattributes( x, {'numeric'}, {'scalar','nonnegative'}, mfilename ) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, S, varargin{:} );

reconDir        = p.Results.recondir;
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

    % Load Sweep Stack
    R(iStk) = load_untouch_nii( S(iStk).rltAbFile(1:end-7), '_swp3d_apod.nii.gz' );
    
    % Load Mask
    N = load_untouch_nii( strcat( S(iStk).maskHeartFile(1:end-7), '_swp3d_apod.nii.gz' ) );
    M{iStk} = single( N.img );
    clear MAT N
    
end


%% Reslice Sweep Volumes into 4-D

for iStk = 1:nStack 
    
    % Reslice Configuration
    nX         = size( R(iStk).img, 1 );
    nY         = size( R(iStk).img, 2 );
    numSwpLoca = max( P(iStk).Sweep.swpWindows(:) ); %TODO: change to size( R(iStk).img, 3 ); so compatible with sweep_window_filter.m ?
    nDyn       = 64;
    nSlices    = numSwpLoca / nDyn;
    
    binWidthSlices = numSwpLoca / nSlices;

    if ~( isreal( binWidthSlices ) && rem( binWidthSlices ,1)==0 )
        error( ['Number of slices does not give integer bin width. Possible bin widths = ' num2str(divisors(numSwpLoca)) ] );
    end
    
    % Init Slice Bins
    edgesSlices = 1:binWidthSlices:numSwpLoca+binWidthSlices;
    binsSlices  = discretize( 1:numSwpLoca, edgesSlices );
    
    % Reslice Sweep Volume
    R_resliced(iStk).img = [];
    
    for iS = 1:nSlices

        currentBinRange = edgesSlices(iS):edgesSlices(iS+1)-1;
        
        R_resliced(iStk).img( :,:,iS,: ) = R(iStk).img( :,:,currentBinRange );

    end
    
end


%% Save NIfTI

for iStk = 1:nStack
    
    % Update img
    R(iStk).img = R_resliced(iStk).img;
    
    % Update Header
    % TODO: do I need to update more fields? Affine?
    nSlices = size( R_resliced(iStk).img, 3 );
    nDyn    = size( R_resliced(iStk).img, 4 );
    R(iStk).hdr.dime.dim([4,5]) = [nSlices, nDyn];     
    
    % Save Resliced Sweep
    S(iStk).rltBinnedAbFile = fullfile( dataDir, sprintf( '%s_rlt_ab_swp_resliced.nii.gz', S(iStk).desc ) );
    save_untouch_nii( R(iStk), S(iStk).rltBinnedAbFile );
    
    % Save Resliced Mask
    % TODO:
    
end


% recon_reslice_sweep(...)
end
