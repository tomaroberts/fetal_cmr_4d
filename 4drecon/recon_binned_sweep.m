function S = recon_binned_sweep( S, varargin )
%RECON_BINNED_SWEEP  discretized binning of sweep data for 4D reconstruction.
%
%   S = RECON_BINNED_SWEEP( reconDir ) uses data structure S from 
%   cardsync_sweep to bin Sweep data into 4D volumes with discrete frames
%   and cardiac phases

%   TAR   (t.roberts@kcl.ac.uk)


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
    
end
    

%% Apodize Sweep Windows
% Filter out artefacted images at start/end of sweep window

apodizationLength = 0; % 0 = no apodization

for iStk = 1:nStack
    
    % TODO: Decide where this filtering comes into the pipeline. Here or in
    % ktrecon code?
    [ R(iStk).img, P(iStk).Sweep ] = sweep_window_filter( R(iStk).img, P(iStk).Sweep, apodizationLength );
    
end


%% Bin Sweep Data 

for iStk = 1:nStack
    
    % Get Cardiac Phases
    thetaFrameSwpBins = [];
    thetaFrameSwpBins = cell2mat( S(iStk).thetaFrameSwpBins );
    thetaFrameSwpBins = thetaFrameSwpBins(:);
    
    
    % Binning Configuration
    nX         = size( R(iStk).img, 1 );
    nY         = size( R(iStk).img, 2 );
    numSwpLoca = max( P(iStk).Sweep.swpWindows(:) ); %TODO: change to size( R(iStk).img, 3 ); so compatible with sweep_window_filter.m ?
    nSlices    = 8;
    nPhases    = 64;

    binWidthSlices = numSwpLoca / nSlices;

    if ~( isreal( binWidthSlices ) && rem( binWidthSlices ,1)==0 )
        error( ['Number of slices does not give integer bin width. Possible bin widths = ' num2str(divisors(numSwpLoca)) ] );
    end
    
    
    % Init Slice and Cardiac Phase Bins
    edgesSlices = 1:binWidthSlices:numSwpLoca+binWidthSlices;
    edgesPhases = linspace(0,2*pi,nPhases+1);

    binsSlices  = discretize( 1:numSwpLoca, edgesSlices );
    binsPhases  = discretize( thetaFrameSwpBins', edgesPhases );
    bins        = [1:numSwpLoca; binsSlices; binsPhases]';

    
    % Bin Sweep Volume into Discrete Slices and Cardiac Phases
    R_binned(iStk).img = zeros( nX, nY, nSlices, nPhases );
    binLogger          = NaN( nSlices ,nPhases );
    
    for iS = 1:nSlices

        for iP = 1:nPhases

            framesInCurrentBin          = find( bins(:,2)==iS & bins(:,3)==iP );
            R_binned(iStk).img( :,:,iS,iP )   = mean( R(iStk).img( :,:,framesInCurrentBin ), 3 ); %TOCHECK: is average correct?
            S(iStk).binLogger( iS, iP ) = numel( framesInCurrentBin );
            clear framesInCurrentBin

        end

    end


    % Save Plot of Slice/Cardiac Phase Bin Distribution
    figure; imagesc( S(iStk).binLogger' );
    axis image; colormap('hot'); colorbar;
    xlabel('Slice Number');
    ylabel('Frame Number');
        
    hFig = gcf;
    hFig.Name = sprintf( 'sweep_bin_dist_nsl%s_nph%s_%s', num2str(nSlices), num2str(nPhases), S(iStk).desc );
    save_figs( pngfigDir, hFig, matfigDir )
    
    if ~isVerbose
        close( hFig )
    end

end


%% Save NIfTI

for iStk = 1:nStack
    
    % Update img
    R(iStk).img = R_binned(iStk).img;
    
    % Update Header
    % TODO: do I need to update more fields? Affine?
    R(iStk).hdr.dime.dim([4,5]) = [nSlices, nPhases];     
    
    % Save
    S(iStk).rltBinnedAbFile = fullfile( ktreconDir, sprintf( '%s_rlt_ab_swpbin.nii.gz', S(iStk).desc ) );
    save_untouch_nii( R(iStk), S(iStk).rltBinnedAbFile );

end


%% Update S structure

save( fullfile( resultsDir, 'results_cardsync_sweep.mat' ), 'S', '-v7.3' );


% RECON_BINNED_SWEEP(...)
end