
%% recon iter_cardsync

% TODO:
% - get specific k-t block
% - get cardiac phases
% - compare / do minimisation


%% Admin

swpDir      = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Sweep';
dataDir     = 'data';
ktreconDir  = 'ktrecon';
cardsyncDir = 'cardsync';
cineVolDir  = 'cine_vol';


%% Run: transform_cine_vol_to_stack.bash
% - transform cine_vol to stack space, e.g.:
% - mirtk transform-image cine_vol_3d.nii.gz cine_vol_3d_transf.nii.gz -target s21_rlt_ab_swp3d_apod.nii.gz -interp NN
% - or: mirtk transform-image cine_vol.nii.gz cine_vol_transf.nii.gz -target s21_rlt_ab_swp3d_apod.nii.gz -interp NN


%% Load Data
cd('C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Sweep\c_fcmr334_3d_iter_recon_cardsync');

% Image Frames
I = load_untouch_nii( 's21_rlt_ab_swp3d_apod.nii.gz' );

% Heart Mask 
M = load_untouch_nii( 's21_mask_heart.nii.gz' );
M.img = single(M.img);

% 4D cine_vol
% X = load_untouch_nii( 'cine_vol_transf.nii.gz' );
X = load_untouch_nii( 'cine_vol_transf_Bspline.nii.gz' );
X.img = single(X.img);

% cine_vol background
B.img = X.img(:,:,:,1) > 0;

% S
iStk = 1;
load( 'results_cardsync_sweep_resliced.mat' );
thetaFrame = S(iStk).thetaFrameSwpLoc;



%% Auto Crop to cine_vol x/y dimensions

for iLoc = 1:size( X.img,3 )

    if ~any( any( X.img(:,:,iLoc,1) ) )
        isEmptySlice(iLoc) = true;
        
        topRow(iLoc)      = nan;
        bottomRow(iLoc)   = nan;
        leftColumn(iLoc)  = nan;
        rightColumn(iLoc) = nan;
        
    elseif any( any( X.img(:,:,iLoc,1) ) )
        isEmptySlice(iLoc) = false;
        
        [ nonZeroRows, nonZeroColumns ] = find( X.img(:,:,iLoc,1) );
        
        topRow(iLoc)      = min( nonZeroRows(:) );
        bottomRow(iLoc)   = max( nonZeroRows(:) );
        leftColumn(iLoc)  = min( nonZeroColumns(:) );
        rightColumn(iLoc) = max( nonZeroColumns(:) );
    
    end

end

cropRegion = [ min(topRow), max(bottomRow), min(leftColumn), max(rightColumn) ];

% Apply Crop
I.img = I.img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : );
M.img = M.img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : );
X.img = X.img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : ,: );
B.img = B.img( cropRegion(1):cropRegion(2), cropRegion(3):cropRegion(4), : );


%% Pad z for consistency with original k-t blocks

I.img = padarray( I.img, [0,0,12], 0, 'both' );
M.img = padarray( M.img, [0,0,12], 0, 'both' );
X.img = padarray( X.img, [0,0,12,0], 0, 'both' );
B.img = padarray( B.img, [0,0,12], 0, 'both' );

thetaFrame = padarray( thetaFrame, 12, NaN, 'both' );


%% thetaFrame
nFrame    = numel(thetaFrame);
nBlocks   = 32;
blockSize = 32;

figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]); 
hold on;

plot( thetaFrame, 'ok', 'MarkerFaceColor','k' );

xticks( 1:blockSize:nFrame );
grid; grid minor;
axis([1 nFrame+10 -0.1 2*pi+0.1]);

title( [ S(iStk).desc ] );
ylabel('Cardiac Phase (theta)');
xlabel('Frame Index');


%% View

% Full FOV
% implay_RR( [I.img.*B.img I.img.*B.img.*M.img X.img(:,:,:,1) X.img(:,:,:,1).*M.img] );

% k-t block

locStart  = 513;
locRange  = locStart:locStart+blockSize-1;

implay_RR( [I.img(:,:,locRange).*B.img(:,:,locRange) ...
            I.img(:,:,locRange).*B.img(:,:,locRange).*M.img(:,:,locRange) ...
            X.img(:,:,locRange,1) ...
            X.img(:,:,locRange,1).*M.img(:,:,locRange)] );
        
        
%% Compare k-t block with corresponding cine_vol

ktBlockStart  = 513;
ktBlockRange  = ktBlockStart:ktBlockStart+blockSize-1;
ktBlockCentre = ktBlockRange( floor(blockSize/2) );

ktBlockthetaFrame = thetaFrame( ktBlockRange );

% cine_vol at centre of k-t block
implay_RR( X.img(:,:,ktBlockCentre,:) );


%% Interp cardiac frames to same number as in k-t block
nX = size(X.img,1); nY = size(X.img,2); nZ = size(X.img,3);
Xinterp.img = zeros( nX, nY, nZ, blockSize );
interpPoints = linspace( 1,25,blockSize );

for xx = 1:nX
    for yy = 1:nY
        for tt = 1:nZ   
            sigt = squeeze( X.img(xx,yy,tt,:) );
            if any(sigt)
                Xinterp.img(xx,yy,tt,:) = interp1( sigt, interpPoints );
            elseif ~any(sigt)
                Xinterp.img(xx,yy,tt,:) = 0;
            end
        end
    end
end

implay_RR( [I.img(:,:,ktBlockRange).*B.img(:,:,ktBlockRange).*M.img(:,:,ktBlockRange) ...
            squeeze(Xinterp.img(:,:,ktBlockCentre,:)).*M.img(:,:,ktBlockRange)] );


%% Compare Signal Intensities

% Threshold
Iblock.img = I.img(:,:,ktBlockRange).*B.img(:,:,ktBlockRange).*M.img(:,:,ktBlockRange);
BP.img     = ~(Iblock.img < prctile( nonzeros( Iblock.img(:) ), 75 ) );
IBP.img    = Iblock.img .* BP.img;

for iF = 1:numel( ktBlockRange )
%     Sig_I(iF) = mean( nonzeros( I.img(:,:,ktBlockRange(iF)).*B.img(:,:,ktBlockRange(iF)).*M.img(:,:,ktBlockRange(iF)) ) );
    Sig_I(iF) = mean( nonzeros( IBP.img(:,:,iF) ) );
end

for iP = 1:size(Xinterp.img,4)
    Sig_X(iP) = mean( nonzeros( Xinterp.img(:,:,ktBlockCentre,iP).*M.img(:,:,ktBlockCentre) ) );
end

% circshift Sig_X so timepoint 1 = ED
[ ~, idxEndDiastoleFrame ] = max( Sig_X );
Sig_X = circshift( Sig_X, -idxEndDiastoleFrame+1 );

% Calculate vector of cardiac phases for cine_vol
thetaFrameCineVol = linspace( 0,2*pi,blockSize ); % Is 0 == 2pi? Should I do blockSize+1?

% figure('units','normalized','outerposition',[0.2 0.2 0.4 0.4]); 
% subplot(1,2,1);
% hold on;
% plot( Sig_I, '-or' );
% xlabel('k-t block frame number');
% ylabel('Signal [a.u.]');
% 
% subplot(1,2,2);
% hold on;
% plot( Sig_X, '-ob' );
% xlabel('Cardiac phase [frame index]');
% ylabel('Signal [a.u.]');


%% Minimise SoS

KTB      = zeros( blockSize, 3 );
KTB(:,1) = 1:blockSize;
KTB(:,2) = Sig_I;
KTB(:,3) = ktBlockthetaFrame;

% Sort by cardiac phase
[KTBsorted(:,3), KTBsorted(:,1)] = sort( KTB(:,3),1 );
KTBsorted(:,2) = Sig_I( KTBsorted(:,1) );

% Perform Minimisation
thetaI = KTBsorted(:,2)';
thetaC = Sig_X;

for iF = 1:numel( ktBlockRange )
    
    thetaIcirc = circshift( thetaI, iF-1 );
    Diff       = thetaIcirc - thetaC;
    ssd(iF)    = sum(Diff(:).^2);
    
end

[ssd_min, ssd_min_idx] = min(ssd);

KTBsorted_shifted = circshift( KTBsorted, ssd_min_idx-1, 1 );

% Assign new theta values to frames
% TODO: check that idxED is the same for every slice in cine_vol
KTBsorted_shifted(:,4) = thetaFrameCineVol;
thetaFrameNew          = sortrows( KTBsorted_shifted );
Sig_I_New              = KTBsorted_shifted(:,2);
thetaFrameNew          = thetaFrameNew(:,4);


% Perturb R-R interval
Perturbations          = 0.03 * 2*pi * ( rand(blockSize,blockSize) );
thetaFrameNewPerturbed = wrapTo2Pi( repmat( thetaFrameNew, 1, blockSize ) .* Perturbations );
thetaI = Sig_I_New';
thetaC = Sig_X;

for iF = 1:numel( ktBlockRange )

	thetaIperturbed = thetaFrameNewPerturbed(:,iF);
    Diff            = thetaIperturbed' - thetaC;
    ssd_perturb(iF) = sum(Diff(:).^2);
    
end

[ssd_min_p, ssd_min_p_idx] = min(ssd_perturb);


%% Plots
figure('units','normalized','outerposition',[0.1 0.2 0.8 0.4]); hold on;

subplot(1,5,1);
hold on;
plot( Sig_I, '-or' );
xlabel('k-t block frame number');
ylabel('Signal [a.u.]');
title(' k-t Block Frames ');

subplot(1,5,2);
plot( KTB(:,3), KTB(:,2), 'xr' );
xlabel('Cardiac phase [radians]');
ylabel('Signal [a.u.]');
title(' k-t Block Before Minimisation ');

subplot(1,5,3);
plot( Sig_X, '-ob' );
xlabel('Cardiac phase [frame index]');
ylabel('Signal [a.u.]');
title(' cine\_vol model ');

subplot(1,5,4);
plot( ssd, '-ko' );
title(' SSD ');
xlabel('Cardiac phase [frame index]');

subplot(1,5,5); hold on;
plot( thetaFrameCineVol, KTBsorted_shifted(:,2)' ./ mean(KTBsorted_shifted(:,2)'), 'rx' );
plot( thetaFrameCineVol, Sig_X ./ mean(Sig_X), '-ob' );
xlabel('New cardiac phase / frame index');
ylabel('Normalised Signals [a.u.]');
title(' k-t Block After Minimisation ');

% subplot(1,6,6);
% hold on;
% plot( Sig_I_New, '-or' );
% xlabel('k-t block frame number');
% ylabel('Signal [a.u.]');
% title(' k-t Block Frames After Minimisation ');


%% k-t Block thetaFrame Plots

figure('units','normalized','outerposition',[0.2 0.2 0.4 0.4]); 

subplot(1,3,1); hold on;
plot( locRange, ktBlockthetaFrame, 'ok', 'MarkerFaceColor','k' );
grid; grid minor;
axis([locRange(1) locRange(end)+1 -0.1 2*pi+0.1]);
title( 'Original Theta' );
ylabel('Cardiac Phase (theta)');
xlabel('Frame Index');

subplot(1,3,2); hold on;
plot( locRange, thetaFrameNew, 'ok', 'MarkerFaceColor','k' );
grid; grid minor;
axis([locRange(1) locRange(end)+1 -0.1 2*pi+0.1]);
title( 'Shifted Theta' );
ylabel('Cardiac Phase (theta)');
xlabel('Frame Index');

subplot(1,3,3); hold on;
plot( locRange, thetaFrameNewPerturbed, 'ok', 'MarkerFaceColor','k' );
grid; grid minor;
axis([locRange(1) locRange(end)+1 -0.1 2*pi+0.1]);
title( 'Shifted Theta + Perturbations' );
ylabel('Cardiac Phase (theta)');
xlabel('Frame Index');
