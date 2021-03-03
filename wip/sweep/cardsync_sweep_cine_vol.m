function S = cardsync_sweep_cine_vol( S, varargin )
%CARDSYNC_SWEEP_CINE_VOL  use cine_vol to perform cardiac synchronisation
%of k-t SWEEP data
%
%   S = CARDSYNC_SWEEP_CINE_VOL( S ) 
%
%   CARDSYNC_SWEEP_CINE_VOL( ..., 'name', value ) specifies optional input argument. 
%   See code for name-value pairs.
%
%   See also .

%   tar (t.roberts@kcl.ac.uk)


%% Optional Input Argument Default Values

default.reconDir        = pwd;
default.resultsDir      = pwd;
default.ktblocksize     = 32;
default.isVerbose       = false;


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

add_param_fn(  p, 'ktblocksize', default.ktBlockSize, ...
        @(x) validateattributes( x, {'numeric'}, {'scalar','nonnegative'}, mfilename ) );

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, S, varargin{:} );

reconDir        = p.Results.recondir;
resultsDir      = p.Results.resultsdir;
ktBlockSize     = p.Results.ktblocksize;
isVerbose       = p.Results.verbose;


%% Setup

dataDir     = fullfile( reconDir, 'data' );
maskDir     = fullfile( reconDir, 'mask' );
ktreconDir  = fullfile( reconDir, 'ktrecon' );
cineVolDir  = fullfile( reconDir, 'cine_vol' );
pngfigDir   = fullfile( resultsDir, 'figs', 'png' );
matfigDir   = fullfile( resultsDir, 'figs', 'fig' );
isVerbose   = false;

nStack      = numel( S );


%% Load Data

for iStk = 1:nStack

    % Load Params
    MAT = matfile( S(iStk).rltParamFile );
    P(iStk) = MAT.PARAM;

    % Load Sweep Data
    R(iStk) = load_untouch_nii( S(iStk).rltAbSwpApodFile );
    
    % Load Masks
    N = load_untouch_nii( S(iStk).maskHeartSwpApodFile );
    M{iStk} = single( N.img );

    clear N MAT
    
end


%% Blood-pool Threshold

for iStk = 1:nStack

    % Apply Mask
    RM{iStk} = R(iStk).img .* M{iStk};
    
    % Blood-pool Threshold
    BP{iStk} = ~(RM{iStk} < prctile( nonzeros( RM{iStk}(:) ), 75 ) );
    
    % Apply Threshold
    RBP{iStk} = R(iStk).img .* BP{iStk};

end




end % cardsync_sweep_cine_vol(...)