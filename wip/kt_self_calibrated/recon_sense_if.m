

%%% README:
% kt Intermediate Frame reconstruction
%
% - This was my initial attempt at implementing intermediate frame
% SENSE reconstruction, ie: take 3 frames and perform SENSE to use as
% training data for full kt reconstruction.
%
% - The aliases were irregular so this method was difficult.
% - In the end, this was superseded by Lucilio's sliding window method,
% which is based on a similar idea.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








function [ xtTrnInt, noiseCov ] = recon_sense_if( ktAcq, csm, psi, varargin )
%RECON_SENSE_IF   create training data using SENSE intermediate frame reconstruction
%
%   ktTrnIf = RECON_SENSE_IF( ktAcq, csm, psi )
%
%   [ xtTrnIf, noiseCov ] = RECON_SENSE_IF( ktAcq, csm, psi )
%
%   RECON_KTSENSE( ..., 'param', val)
%
%   output:
%       xtTrnIf             reconstructed training data         (complex,   x-y-time)
%       noiseCov            noise covariance matrix             (complex,   channel-channel)
%
%   input: 
%       ktAcq               k-t undersampled k-space data       (complex,   x-y-time-channel)
%       csm                 coil/channel sensitivity maps       (complex,   x-y-1-channel)
%       psi                 noise covariance matrix;            (complex,   channel-channel)
%
%   option parameter-value pairs:
%       Reff                Intermediate frames SENSE factor    (real/positive, scalar)
%       lambda              SENSE Tikhonov regularisation       (real/positive, scalar)
%       regmethod           SENSE regularisation method         (string)
%       TO USE?: lambdaroi           regularisation inside ROI          (real/positive, scalar)
%       removeoversampling  remove frequency oversampling       (logical,   scalar)
%       fn_rmv_os           remove oversampling function        (anonymous function)
%       verbose             output intermediate steps           (logical,   scalar)
%
%   see also: recon_ktsense, recon_xtsense, recon_xfsense

% tar (t.roberts@kcl.ac.uk) // thanks to jfpva (joshua.vanamerom@kcl.ac.uk) 

%% Notes

% Definitely not good yet...
%
% TODO: make sliding window priors, which are centred around intermediate frames
% TODO: SENSE reconstruction currently based on lots of for loops - make more efficient
%


%% Data Dimensions

dimX = 1;  nX = size( ktAcq, dimX );    % frequency-encode
dimY = 2;  nY = size( ktAcq, dimY );    % phase-encode
dimT = 3;  nT = size( ktAcq, dimT );    % time
dimC = 4;  nC = size( ktAcq, dimC );    % channel
dimF = dimT;  nF = nT;


%% Parse Input

default.Reff                = 3;
default.lambda              = 0.1;
default.regmethod           = 'tikhonov';
default.noiseCov            = [];
default.removeoversampling  = false;
default.fn_rmv_os           = @( x ) x((round(size(x,dimX)/4)+1):(size(x,dimX)-round(size(x,dimX)/4)),:,:,:);
default.isVerbose           = false;

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(  p, 'ktAcq', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nX,nY,nT,nC]}, mfilename ) );

addRequired(  p, 'csm', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[NaN,nY,1,nC]}, mfilename ) );
    
addRequired(  p, 'psi', @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nC,nC]}, mfilename ) );

add_param_fn( p, 'Reff', default.Reff, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );
    
add_param_fn( p, 'lambda', default.lambda, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'positive','scalar'}, mfilename ) );
    
add_param_fn( p, 'regmethod', default.regmethod, ...
        @(x) validateattributes( x, {'char'}, ...
        {'vector'}, mfilename ) );

add_param_fn( p, 'noisecov', default.noiseCov, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'size',[nC,nC]}, mfilename ) );
    
add_param_fn( p, 'removeoversampling', default.removeoversampling, ...
        @(x) validateattributes( x, {'logical'}, ...
        {'scalar'}, mfilename ) );

add_param_fn( p, 'fn_rmv_os', default.fn_rmv_os, ...
        @(x) validateattributes( x, {'function_handle'}, ...
        {}, mfilename ) );

add_param_fn( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, ktAcq, csm, psi, varargin{:} );

Reff                = p.Results.Reff;
lambda              = p.Results.lambda;
regmethod           = p.Results.regmethod;
noiseCov            = p.Results.noisecov;
isRmvOversampling   = p.Results.removeoversampling;
fn_rmv_os           = p.Results.fn_rmv_os;
isVerbose           = p.Results.verbose;

assert( ~isreal( ktAcq ), 'ktAcq not complex-valued' )
assert( ~isreal( csm ), 'csm not complex-valued' )
assert( ~isreal( psi ), 'psi not complex-valued' )


%% Setup

if ( isVerbose )
    fprintf( '\n%s()\n\n', mfilename );
end


%% Anon Fns

kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, dimX ), dimY ) );

phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [+1;-1], size(k,1)/2, 1 ) ) ) );


%% Phase Correction

ktAcq = phase_correct( ktAcq );


%% Transform to x-t Space

xtAcq = kt2xt( ktAcq );


%% View Intermediate k-space frames

if ( isVerbose )
    
    Reff = 3; % effective SENSE factor
    coilToView = 1;
    
    figure;
    for ii = 1:Reff
        subplot(1,3,ii);
        frameLines = ii*Reff-Reff+1:ii*Reff;
        imagesc(abs(sum(ktAcq(:,:,frameLines,coilToView),3)),[0,100]);
        colormap(gray), colorbar('off'), axis image;
        title(['Intermediate Frame ' num2str(ii)]);
    end

    % View aliased Intermediate frames
    figure;
    for ii = 1:Reff
        subplot(1,3,ii);
        frameLines = ii*Reff-Reff+1:ii*Reff;
        imagesc(abs(sum(xtAcq(:,:,frameLines,coilToView),3)),[0,0.1]);
        colormap(gray), colorbar('off'), axis image;
        title(['Intermediate Frame ' num2str(ii)]);
    end
    
end


%% Create Intermediate Frames

if rem(nT/Reff,1) ~= 0
    error('Cannot construct intermediate frames. nY not divisible by Reff');
end

ktInt = squeeze( sum( reshape( ktAcq, nX, nY, Reff, nT/Reff, nC ), 3 ) );

xtInt = kt2xt( ktInt );

nI = size(xtInt,3);


%% Build Priors for Tikhonov Regularisation

ktPri = squeeze( sum( reshape( ktAcq, nX, nY, 8, [], nC ), 3 ) );
%TODO: sliding window, so ktPri is same size as ktAcq?

xtPri = kt2xt( ktPri );

% Sum of Squares Priors
xtPri = sqrt( sum( abs(xtPri).^2 , dimC ) );

% Duplicate to match size of Intermediate Frames
%TODO: make sliding window priors
nP = size(xtPri,3);
for pp = 1:nP
    temp(:,:,1+Reff*(pp-1):Reff*pp) = repmat(xtPri(:,:,pp),[1,1,3]);
end
xtPri = temp; clear temp;


%% Trim Data FOVs
% nY must be divisible by Reff

nY_mod = mod(nY,Reff);
nY_aliased = (nY-nY_mod) / Reff;

xtIntRed = xtInt(:,1:nY_aliased,:,:);
xtPriRed = xtPri(:,1:nY-nY_mod,:);
csmRed   = csm(:,1:nY-nY_mod,:,:);


%% Perform SENSE reconstruction using intermediate frames

% regmethod = 'off';
regmethod = 'tikhonov';

tol = 1e-9; % pinv tolerance
inv_psi = inv(psi);

xtTrnInt = zeros( nX, nY, nI );

switch regmethod
    
    % No Regularisation
    case { 'off' , 'none' }

        for ii = 1:nI
            
            for y = 1:nY_aliased
                
                y_idx = y:nY_aliased:nY_aliased*Reff;
                
                for x = 1:nX
                    
                    % Sensitivities for current voxels
                    S = transpose(reshape(csmRed(x,y_idx,1,:),Reff,[]));
                    
                    % SENSE Encoding Matrix
                    unmix = pinv( S'*inv_psi*S , tol ) * S' * inv_psi;
                    
                    % Unliased Images
                    xtTrnInt(x,y_idx,ii) = unmix * reshape(xtIntRed(x,y,ii,:),[],1);
                    
                    
                end
            end
            
            disp(['Finished intermediate frame number: ' num2str(ii)]);
            
        end
    
    % Tikhonov Regularisation
    case 'tikhonov'
        
        A = eye(Reff);
        AinvA = A'*A;
        
        for ll = 1:numel(lambda)
        
            for ii = 1:nI
            
                for y = 1:nY_aliased

                    y_idx = y:nY_aliased:nY_aliased*Reff;

                    for x = 1:nX

                        % Sensitivities for current voxels
                        S = transpose(reshape(csmRed(x,y_idx,1,:),Reff,[]));

                        % SENSE Encoding Matrix
                        unmix = pinv( S'*inv_psi*S + lambda(ll).*AinvA  , tol ) * S' * inv_psi;

                        % Difference between Unaliased and Prior Images
                        diffImPri = reshape(xtIntRed(x,y,ii,:),[],1) - S*transpose(xtPriRed(x,y_idx,ii));
                        
                        % Unliased Images
                        xtTrnInt(x,y_idx,ii,ll) = transpose(xtPriRed(x,y_idx,ii)) + unmix * diffImPri;               

                    end
                    
                end
            
            disp(['Finished intermediate frame number: ' num2str(ii)]);
                
            end
        
        if numel(lambda) > 1
            disp(['Finished reconstruction number: ' num2str(ll) ' ... lambda = ' num2str(lambda(ll))]);
        end
            
        end
        
end



end  % recon_sense_if(...)