% Run: C:\Users\tr17\Documents\CODE_Projects\fetal_cmr_4d\wip\sweep\kt_block_recon_iter_cardsync.m

% Then:

%% Circshift - align frames with cine_vol model
KTB      = zeros( blockSize, 3 );
KTB(:,1) = 1:blockSize;             % thetaFrame indices
KTB(:,2) = Sig_I;                   % frame Signal intensities
KTB(:,3) = ktBlockthetaFrame;       % thetaFrame

% Sort by cardiac phase
KTB = sortrows( KTB, 3 );
Sig_I_sorted = KTB(:,2)';
ssd = [];

% Perform Minimisation
for iF = 1:numel( ktBlockRange )
    
    Sig_I_circ = circshift( Sig_I_sorted, iF-1 );
    Diff       = Sig_I_circ - Sig_X;
    ssd(iF)    = sum(Diff(:).^2);
    
end

[ssd_min, ssd_min_idx] = min( ssd );

% Shift thetaFrame
thetaShift = 2*pi * ((ssd_min_idx-1)/blockSize); % -1 correct?
KTB        = sortrows( KTB, 1 ); % resort by frame order
KTB(:,3)   = wrapTo2Pi( KTB(:,3) + thetaShift );

ktBlockthetaFrameShifted = KTB(:,3);


%% Perturb R-R interval

% Setup interpolated cine_vol signal curve
ntheta       = 1000;
theta        = linspace(0,2*pi,ntheta);
Sig_X_interp = interp1( thetaFrameCineVol, Sig_X, theta );

% Set R-R interval Perturbations
deltaTheta    = mean( wrapTo2Pi( diff( ktBlockthetaFrame ) ) );
Perturbations = deltaTheta .* [0.5 0.95 0.97 0.98 0.99 1 1.01 1.02 1.03 1.05 1.5];
% Perturbations = deltaTheta .* [0.5:0.05:1.5];
% Perturbations = deltaTheta .* [0.98:0.005:1.02];

% Extend KTB
KTB(:,1:2,1:numel(Perturbations)) = repmat( KTB(:,1:2), [1,1,numel( Perturbations)] );

% Apply Perturbations to thetaFrame
for iP = 1:numel( Perturbations )
    KTB(:,3,iP) = ktBlockthetaFrameShifted(1) : Perturbations(iP) : ktBlockthetaFrameShifted(1)+(Perturbations(iP)*(blockSize-1));
end
KTB(:,3,:) = wrapTo2Pi( KTB(:,3,:) );

% Plot
figure;
for iP = 1:numel( Perturbations )
    subplot(4,3,iP); hold on;
    plot( locRange, KTB(:,3,iP), 'ok', 'MarkerFaceColor','k' );
    grid; grid minor;
    axis([locRange(1) locRange(end)+1 -0.1 2*pi+0.1]);
    title(['Perturbation = ' num2str(deltaTheta - Perturbations(iP))]);
end


% Minimisation
figure;
ssd_pert = [];

for iP = 1:numel( Perturbations )
% for iP = 5

    KTBtemp      = sortrows( KTB(:,:,iP), 3 );
    Sig_I_sorted = KTBtemp(:,2)';
    thetaI       = KTBtemp(:,3)';
    
    Sig_I_sorted = Sig_I_sorted ./ mean( Sig_I_sorted );
    
    % Find interpolated indices of thetaI
    idxThetaInterp(iP,:) = dsearchn( theta', thetaI' );
    
    % cine_vol interpolated signal
%     Sig_X_Matrix(iP,:) = Sig_X_interp( idxThetaInterp(iP,:) );
    Sig_X_Matrix(iP,:) = Sig_X_interp( idxThetaInterp(iP,:) ) ./ mean( Sig_X_interp( idxThetaInterp(iP,:) ) );

    
    % SSD
    Diff         = Sig_I_sorted - Sig_X_Matrix(iP,:);
    ssd_pert(iP) = sum(Diff(:).^2);
   
    % Residuals
    subplot(4,3,iP); hold on;
    plot( thetaI,  Diff, 'ko' );
%     plot( thetaI, Sig_I_sorted ./ mean(Sig_I_sorted), 'rx' );
    title(['Perturbation = ' num2str(deltaTheta - Perturbations(iP))]);
    
    
%     % Plot
%     subplot(4,3,iP); hold on;
% %     plot( theta,  Sig_X_interp ./ mean(Sig_X_interp), 'b' );
% %     plot( thetaI, Sig_I_sorted ./ mean(Sig_I_sorted), 'rx' );
%     plot( thetaI,  Sig_X_Matrix(iP,:), 'b' );
%     plot( thetaI, Sig_I_sorted, 'rx' );
%     title(['Perturbation = ' num2str(deltaTheta - Perturbations(iP))]);
    
end

[ssd_pert_min, ssd_pert_min_idx] = min( ssd_pert );

% Shift thetaFrame
KTBtemp = sortrows( KTB(:,:,ssd_pert_min_idx), 1 ); % resort by frame order
% KTBtemp = sortrows( KTB(:,:,5), 1 ); % resort by frame order
ktBlockthetaFramePerturbed = KTBtemp(:,3);


%% Plot SSDs
figure; hold on; 
subplot(1,2,1); 
plot( ssd, 'r' );
title('Frame to Cine Alignment');
subplot(1,2,2); 
plot( Perturbations, ssd_pert, '-bo' );
title('R-R Perturbations');


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
plot( locRange, ktBlockthetaFrameShifted, 'ok', 'MarkerFaceColor','k' );
grid; grid minor;
axis([locRange(1) locRange(end)+1 -0.1 2*pi+0.1]);
title( 'Shifted Theta' );
ylabel('Cardiac Phase (theta)');
xlabel('Frame Index');

subplot(1,3,3); hold on;
plot( locRange, ktBlockthetaFrameShifted, 'ok', 'MarkerSize', 9);
plot( locRange, ktBlockthetaFramePerturbed, 'ok', 'MarkerFaceColor','k' );
grid; grid minor;
axis([locRange(1) locRange(end)+1 -0.1 2*pi+0.1]);
title( 'Shifted Theta + Perturbations' );
ylabel('Cardiac Phase (theta)');
xlabel('Frame Index');
