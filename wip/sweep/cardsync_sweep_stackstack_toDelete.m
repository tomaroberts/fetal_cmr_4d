
function S = cardsync_sweep_interstack( S, varargin )
%CARDSYNC_SWEEP_INTRASTACK  Align cardiac phases of Sweep volumes.
%
%   S = CARDSYNC_SWEEP_INTERSTACK( S ) 
%
%   CARDSYNC_SWEEP_INTERSTACK( ..., 'name', value ) specifies optional input argument. 
%   See code for name-value pairs.
%
%   See also PREPROC.

%   tar (t.roberts@kcl.ac.uk)


%% Optional Input Argument Default Values

default.reconDir        = pwd;
default.resultsDir      = pwd;       % path of directory to save results
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

add_param_fn(   p, 'verbose', default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, S, varargin{:} );

reconDir        = p.Results.recondir;
resultsDir      = p.Results.resultsdir;
isVerbose       = p.Results.verbose;


%% Setup

dataDir     = fullfile( reconDir, 'data' );
maskDir     = fullfile( reconDir, 'mask' );
ktreconDir  = fullfile( reconDir, 'ktrecon' );
pngfigDir   = fullfile( resultsDir, 'cardsync', 'figs', 'png' );
matfigDir   = fullfile( resultsDir, 'cardsync', 'figs', 'fig' );
isVerbose   = false;

nStack      = numel( S );


%% Load Windowed Sweep Data

for iStk = 1:nStack

    % Load Params
    M = matfile( S(iStk).rltParamFile );
    P(iStk) = M.PARAM;

    % Load NIfTI
    R(iStk) = load_untouch_nii( S(iStk).rltAbFile );
    
    clear M
    
end
    

%% Apodize Sweep Windows
% Filter out artefacted images at start/end of sweep window

apodizationLength = 0; % 0 = no apodization

for iStk = 1:nStack
    
    nX   = size( R(iStk).img, 1 );
    nY   = size( R(iStk).img, 2 );
    nZ   = size( R(iStk).img, 3 );
    nT   = size( R(iStk).img, 4 );
    
    % TODO: Decide where this filtering comes into the pipeline. Here or in
    % ktrecon code?
    [ R(iStk).img, P(iStk).Sweep ] = sweep_window_filter( R(iStk).img, P(iStk).Sweep, apodizationLength );
    
    % Adjust 4d parameters to 3d parameters
    R(iStk).hdr.dime.dim(1) = 3;
    R(iStk).hdr.dime.dim(4) = size( R(iStk).img,3 );
    R(iStk).hdr.dime.dim(5) = 1;

    % nii3d.hdr.dime.pixdim(1) = 0; % Not sure necessary? 0 or 1?
    R(iStk).hdr.dime.pixdim(4) = R(iStk).hdr.dime.pixdim(4) / nT;

    % Save
    save_untouch_nii( R(iStk), [S(iStk).rltAbFile(1:end-7) '_swp3d_apod.nii.gz' ] );    
    
end


%% Stack-stack cardsync
% Maximise cross-correlation


%% Load Masks

% for iStk = 1:nStack
for iStk = 1
   
    N = load_untouch_nii( strcat( S(iStk).maskHeartFile(1:end-7), '_swp3d_apod.nii.gz' ) );
    M{iStk} = single( N.img );

    clear N
    
end


%% Blood-pool Threshold

RM = R(iStk).img .* M{iStk};

for iStk = 1
   
    BP{iStk} = ~(RM < prctile( nonzeros(RM(:)), 75 ) ); % percentile method
    
end

% iS = 530;
% imtar(RM(:,:,iS));
% imtar(BP{1}(:,:,iS));

M = BP;


%% Calculate Mean Intensity in Heart Region

% for iStk = 1:nStack
for iStk = 1
    
    numSwpLoca = max( P(iStk).Sweep.swpWindows(:) );
    
    for iLoc = 1:numSwpLoca     
        
        if ~any( M{iStk}(:,:,iLoc), 'all' )
            S(iStk).meanSignalHeart(iLoc) = NaN;
        else
            S(iStk).meanSignalHeart(iLoc) = mean( nonzeros( R(iStk).img(:,:,iLoc) .* M{iStk}(:,:,iLoc) ) );
        end
                 
    end    

end

figure; 
plot(S(iStk).meanSignalHeart,'.-k');
AX = axis; axis([1 numSwpLoca AX(3) AX(4)]);


%% Sort Cardiac Phases

% % for iStk = 1:nStack
% for iStk = 1
%     
%     thetaFrameSwpBins = [];
%     thetaFrameSwpBins = cell2mat( S(iStk).thetaFrameSwpBins );
%     thetaFrameSwpBins = thetaFrameSwpBins(:);
%    
%     [thetaFrameSorted, itF ] = sort( thetaFrameSwpBins );
%     
%     figure;
%     plot( S(iStk).meanSignalHeart( itF ) );
%     
% end


%% Bin Cardiac Phases to Estimate Cardiac Volume Curve


% for iStk = 1:nStack
for iStk = 1

    thetaFrameSwpBins = [];
    thetaFrameSwpBins = cell2mat( S(iStk).thetaFrameSwpBins );
    thetaFrameSwpBins = thetaFrameSwpBins(:);
    
    % Binning Configuration
    nX         = size( R(iStk).img, 1 );
    nY         = size( R(iStk).img, 2 );
    numSwpLoca = max( P(iStk).Sweep.swpWindows(:) ); %TODO: change to size( R(iStk).img, 3 ); so compatible with sweep_window_filter.m ?
    nPhases    = 32;   
    
    % Init Slice and Cardiac Phase Bins
    edgesPhases = linspace(0,2*pi,nPhases+1);
    binsPhases  = discretize( thetaFrameSwpBins', edgesPhases );

    % Mean Signal over Frames in Each Phase Range
    meanSignalIntensity = [];
    
    for iP = 1:nPhases
        
        framesInCurrentPhaseRange = find( binsPhases == iP );
        
        meanSignalIntensity(iP) = mean( nonzeros( R(iStk).img(:,:,framesInCurrentPhaseRange) .* ...
                                                  M{iStk}(:,:,framesInCurrentPhaseRange) ) );
                                            
        clear framesInCurrentPhaseRange
        
    end

end

figure; 
plot(edgesPhases(1:end-1), meanSignalIntensity,'.-k');
xlabel('Cardiac Phase [radians]');
ylabel('Mean Signal Intensity');



% cardsync_sweep_interstack(...)
end
