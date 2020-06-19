%% Voxel-wise comparison between PC-SPGR and PC-bSSFP
%
% - script generates histograms or AP/FH/LR/magnitude velocity values
% within the pipes of the flow phantom
% - FlowPhantom9 scans used for analysis - these were matched resolutions
% and FOVs
%
%


%% Data Choices
qflowDataDir = 'data_raw_rigid_to_ap';
bssfpDataDir = 'data_rigid_to_qflow';

numBssfpStacks = 5;
extVelDir      = '_noSag';
useIsotropic   = true;  % use isotropic volumes
% useIsotropic   = false;
useInterp      = false; % use resampled, isotropic volumes
useErodedMask  = false;


%% Admin
fp9Dir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9';

qflowReconDir = [fp9Dir '\QFlow'];
bssfpReconDir = [fp9Dir '\PC-bSSFP'];


%% Load QFlow data
cd(qflowReconDir);
cd(qflowDataDir);

switch qflowDataDir
        
        % QFlow from .nii
        case 'data'

            velFiles = dir('*_vel.nii.gz');
            % velFiles = dir('*_vel_*scaled.nii.gz');
            velFiles = {velFiles([1,2,3]).name};
            for ii = 1:numel(velFiles)    
                nii = load_nii(velFiles{ii});   
                VEL(:,:,:,ii) = nii.img;
            end

        % QFlow from .nii, with FFD registration
        case 'data_ffd'
            
            velFiles = dir('*_vel_*scaled.nii.gz');
            velFiles = {velFiles([1,2,3]).name};
            for ii = 1:numel(velFiles)    
                nii = load_nii(velFiles{ii});   
                VEL(:,:,:,ii) = nii.img;
            end
            
        % QFlow from .raw, with no registration
        case 'data_raw'

            phFiles = dir('*_ph_corr.nii.gz');
            phFiles = {phFiles([1,2,3]).name};
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);
            
        % QFlow from .raw, with rigid registration
        case 'data_raw_rigid'
            
            phFiles = dir('*_ph_corr_rigid.nii.gz');
            phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);
            
        % QFlow from .raw, with FFD registration
        case 'data_raw_ffd'

            phFiles = dir('*_ph_corr_ffd.nii.gz');
            phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);
            
        % QFlow from .raw, rigid registered to QFlow AP stack
        case 'data_raw_rigid_to_ap'
            
            if ( useIsotropic && useInterp )
                cd isotropic_interp_0p5
                phFiles = dir('*_ph_corr_rigid_isotropic.nii.gz');
            elseif useIsotropic
                cd isotropic
                phFiles = dir('*_ph_corr_rigid_isotropic.nii.gz');
            else
                phFiles = dir('*_ph_corr_rigid.nii.gz');
            end
            
            phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);
            
        % QFlow from .raw, affine registered to QFlow AP stack
        case 'data_raw_affine_to_ap'
        
            phFiles = dir('*_ph_corr_affine.nii.gz');
            phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);       
        
end

vq = double(VEL);



%% Load PC-bSSFP data
cd(bssfpReconDir);
cd(bssfpDataDir);
cd([num2str(numBssfpStacks) '_phase_stacks_recon' extVelDir]);

switch bssfpDataDir
        
        % PC-bSSFP transformed to axial stack
        case 'data_trans'
            
            nii = load_nii('ap_vel_trans.nii.gz');
            vb(:,:,:,1) = double(nii.img);
            nii = load_nii('fh_vel_trans.nii.gz');
            vb(:,:,:,2) = -1.*double(nii.img);
            nii = load_nii('rl_vel_trans.nii.gz');
            vb(:,:,:,3) = double(nii.img);
            
            
        % PC-bSSFP rigid registered to axial bSSFP stack
        case 'data_rigid'
            
            nii = load_nii('ap_vel_rigid.nii.gz');
            vb(:,:,:,1) = double(nii.img);
            nii = load_nii('fh_vel_rigid.nii.gz');
            vb(:,:,:,2) = -1.*double(nii.img);
            nii = load_nii('rl_vel_rigid.nii.gz');
            vb(:,:,:,3) = double(nii.img);
            
            
        % PC-bSSFP FFD registered to axial bSSFP stack
        case 'data_ffd'
            
            nii = load_nii('ap_vel_ffd.nii.gz');
            vb(:,:,:,1) = double(nii.img);
            nii = load_nii('fh_vel_ffd.nii.gz');
            vb(:,:,:,2) = -1.*double(nii.img);
            nii = load_nii('rl_vel_ffd.nii.gz');
            vb(:,:,:,3) = double(nii.img);
            
            
        % PC-bSSFP rigid registered to axial QFlow stack
        case 'data_rigid_to_qflow'
            
            if ( useIsotropic && useInterp )
                cd isotropic_interp_0p5
                velFiles = dir('*_rigid_isotropic.nii.gz');                
                velFiles = {velFiles([1,2,3]).name}; %re-order: ap/fh/rl
                for ii = 1:numel(velFiles)    
                    nii = load_nii(velFiles{ii});   
                    vb(:,:,:,ii) = double(nii.img);
                end
            elseif useIsotropic
                cd isotropic
                velFiles = dir('*_rigid_isotropic.nii.gz');                
                velFiles = {velFiles([1,2,3]).name}; %re-order: ap/fh/rl
                for ii = 1:numel(velFiles)    
                    nii = load_nii(velFiles{ii});   
                    vb(:,:,:,ii) = double(nii.img);
                end
            else
                nii = load_nii('ap_vel_rigid.nii.gz');
                vb(:,:,:,1) = double(nii.img);
                nii = load_nii('fh_vel_rigid.nii.gz');
                vb(:,:,:,2) = -1.*double(nii.img);
                nii = load_nii('rl_vel_rigid.nii.gz');
                vb(:,:,:,3) = double(nii.img);
            end
            



        % PC-bSSFP affine registered to axial QFlow stack
        case 'data_affine_to_qflow'

            nii = load_nii('ap_vel_affine.nii.gz');
            vb(:,:,:,1) = double(nii.img);
            nii = load_nii('fh_vel_affine.nii.gz');
            vb(:,:,:,2) = -1.*double(nii.img);
            nii = load_nii('rl_vel_affine.nii.gz');
            vb(:,:,:,3) = double(nii.img);
            
            
end


%% Create Velocity Magnitude Volumes
vq3d = sqrt(vq(:,:,:,1).^2 + vq(:,:,:,2).^2 + vq(:,:,:,3).^2);
vb3d = sqrt(vb(:,:,:,1).^2 + vb(:,:,:,2).^2 + vb(:,:,:,3).^2);


%% Load Water Pipe Mask
cd(qflowReconDir);
cd mask

if ( useIsotropic && useInterp )
    P = load_nii('mask/inside_pipes_mask3_isotropic_interp_0p5.nii.gz');
elseif useIsotropic
    P = load_nii('mask/inside_pipes_mask3_isotropic.nii.gz');
else
    P = load_nii('mask/inside_pipes_mask3.nii.gz');
    % S = load_nii('spherical_part_mask.nii.gz');
end

M = double(P.img);

% Erode mask
if useErodedMask == true
    SE = strel('sphere',1); Me = imerode(M,SE);
    M = Me;
end

%% View Velocity Images
implay_RR([M.*vq(:,:,:,1), M.*vq(:,:,:,2), M.*vq(:,:,:,3), M.*vq3d;...
           M.*vb(:,:,:,1), M.*vb(:,:,:,2), M.*vb(:,:,:,3), M.*vb3d]...
           ,'jet',[-10,10]);


%% Flow in Pipe ROIs
% drawn on ../data_raw_rigid_to_ap/ap_mag_rigid.nii.gz
roiDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow\roi';
cd(roiDir);

if ( useIsotropic && useInterp )
    cd isotropic_interp_0p5
    roiFiles = dir('roi*.nii.gz');
elseif useIsotropic
    cd isotropic
    roiFiles = dir('roi*.nii.gz');
else
    roiFiles = dir('roi*.nii.gz');
end

ROI = zeros(size(vq,1),size(vq,2),size(vq,3),numel(roiFiles));

for ii = 1:numel(roiFiles)
    nii = load_nii(roiFiles(ii).name);
    
    ROI(:,:,:,ii) = M .* double(nii.img); % NB: apply threshold mask to ROIs to remove background/pipe voxels
%     ROI(:,:,:,ii) = double(nii.img); % don't apply threshold mask
    
    numVoxInROI(ii)   = numel(nonzeros( ROI(:,:,:,ii) ));
    iROI{ii}          = find( ROI(:,:,:,ii) );
    
    vq3d_ROI{ii}      = vq3d( iROI{ii} );
    vb3d_ROI{ii}      = vb3d( iROI{ii} );
    
    vq3d_ROI_mean(ii) = mean( vq3d_ROI{ii} );
    vb3d_ROI_mean(ii) = mean( vb3d_ROI{ii} );
end


% convert current velocity through ROI to flow
voxRes = nii.hdr.dime.pixdim(2:4);
voxArea = voxRes(1) * voxRes(2); % mm^2

vel2flow = @( v, a ) ( v ) .* ( a * 1e-2 );  % ml/s

for ii = 1:numel(roiFiles)

    roiArea = voxArea * numel(iROI{ii});
    
    velMean = vq3d_ROI_mean(ii);
    fq3d_ROI_mean(ii,1) = vel2flow( velMean, roiArea ); % [ml/s]
    
    velMean = vb3d_ROI_mean(ii);
    fb3d_ROI_mean(ii,1) = vel2flow( velMean, roiArea ); % [ml/s]

end

% Benchtest
MeanFlow.BenchTest = (0.555 + 0.932) / 2; % [mL/s] - Mean of two bench tests

% Correlation plot
figure; hold on;
plot(fq3d_ROI_mean(1:end),fb3d_ROI_mean(1:end),'.k','Markersize',20);
plot(MeanFlow.BenchTest,MeanFlow.BenchTest,'r.','Markersize',20);
plot(0:0.1:2,0:0.1:2,'k--');
xlabel('QFlow [ml/s]');
ylabel('bSSFP [ml/s]');
axis([0.3 1.1 0.3 1.1]);
legend('Flow MRI','Flow Bench Measurement','Unity','Location','SouthEast');

% Histograms
% figure('units','normalized','outerposition',[0 0 1 1]);
% ctr = 1;
% for ii = 1:10
%     subplot(2,10,ii);
%     histogram(vq3d_ROI{ii},20); axis([0 10 0 30]); title(['QFLOW ROI no. ' num2str(ctr)]);
%     subplot(2,10,ii+10);
%     histogram(vb3d_ROI{ii},20); axis([0 10 0 30]); title(['bSSFP ROI no. ' num2str(ctr)]);
%     ctr = ctr+1;
% end

% Save for GraphPad

MeanFlow.SPGR  = mean(fq3d_ROI_mean);
MeanFlow.bSSFP = mean(fb3d_ROI_mean);
StdFlow.SPGR   = std(fq3d_ROI_mean);
StdFlow.bSSFP  = std(fb3d_ROI_mean);

save('MeanFlows_20_ROIs.mat','fq3d_ROI_mean','fb3d_ROI_mean', ...
     'MeanFlow', 'StdFlow');


%% View Flow Images
% implay_RR([M.*fq3d;...
%            M.*fb3d]...
%            ,'jet',[-.01,.01]);
       

%% Plot Histograms

lowerThresh = -20;
upperThresh = 20;
binSize = 0.2;

bkrdClip = 0.025;

for ii = 1:3
    M_vq(:,:,:,ii) = M.*vq(:,:,:,ii);
    M_vb(:,:,:,ii) = M.*vb(:,:,:,ii);
end

M_vq3d = M.*vq3d;
M_vb3d = M.*vb3d;

% Set background/really low values to NaN
% - NaNs excluded from histogram function
% - use QFlow velocity magnitude for thresholding
nanMask = double(M_vq3d>bkrdClip); nanMask(nanMask<1)=NaN;

for ii = 1:3
    tempVect = M_vq(:,:,:,ii).*nanMask; QVhist(:,ii) = tempVect(:);
    tempVect = M_vb(:,:,:,ii).*nanMask; BVhist(:,ii) = tempVect(:);
end

tempVect = M_vq3d .* nanMask; QVhist(:,4) = tempVect(:);
tempVect = M_vb3d .* nanMask; BVhist(:,4) = tempVect(:);


%%% Components
figure;
for ii = 1:3
    
    subplot(2,2,ii);
    histogram(QVhist(:,ii),'binEdges',lowerThresh:binSize:upperThresh,'Normalization','Probability');
    hold on;
    histogram(BVhist(:,ii),'binEdges',lowerThresh:binSize:upperThresh,'Normalization','Probability');
    
    legend('PC-SPGR','PC-bSSFP');
    xlabel('Velocity [cm/s]');
    ylabel('Counts');
    

    if ii == 1
        title('AP');
    elseif ii == 2
        title('FH');
    elseif ii == 3
        title('RL');
    end
    
    % For Bland-Altmans
    QnonNaNidx = find(~isnan(QVhist(:,ii)));
    BnonNaNidx = find(~isnan(BVhist(:,ii)));
    
    velsBlandAltman{ii}(:,1) = QVhist(QnonNaNidx,ii);
    velsBlandAltman{ii}(:,2) = BVhist(BnonNaNidx,ii);
    
end


%%% Magnitude
subplot(2,2,4);
histogram(QVhist(:,4),'binEdges',0:binSize/2:upperThresh,'Normalization','Probability');
hold on;
histogram(BVhist(:,4),'binEdges',0:binSize/2:upperThresh,'Normalization','Probability');

legend('PC-SPGR','PC-bSSFP');
xlabel('Velocity [cm/s]');
ylabel('Counts');
title('Magnitude');

% For Bland-Altmans
QnonNaNidx = find(~isnan(QVhist(:,4)));
BnonNaNidx = find(~isnan(BVhist(:,4)));

velsBlandAltman{4}(:,1) = QVhist(QnonNaNidx,4);
velsBlandAltman{4}(:,2) = BVhist(BnonNaNidx,4);


%% Old stuff for PDFs
% figure; hold on;
% Qpdca = fitdist(QVhist(QVhist<10),'Normal');
% Bpdca = fitdist(BVhist(BVhist<10),'Normal');

% Qpdca = fitdist(QVhist,'Normal');
% Bpdca = fitdist(BVhist,'Normal');
% 
% figure; hold on;
% plot([0:0.25:10],pdf(Qpdca,[0:0.25:10]));
% plot([0:0.25:10],pdf(Bpdca,[0:0.25:10]));
% legend('PC-SPGR','PC-bSSFP');
% xlabel('Velocity [cm/s]');


% %% For Bland-Altmans
% - For some reason this code broke...
% BA.ap(:,1)  = QVhist(~isnan(BVhist(:,1))); %ap
% BA.ap(:,2)  = BVhist(~isnan(BVhist(:,1)));
% BA.fh(:,1)  = QVhist(~isnan(BVhist(:,2))); %fh
% BA.fh(:,2)  = BVhist(~isnan(BVhist(:,2)));
% BA.rl(:,1)  = QVhist(~isnan(BVhist(:,3))); %rl
% BA.rl(:,2)  = BVhist(~isnan(BVhist(:,3)));
% BA.mag(:,1) = QVhist(~isnan(BVhist(:,4))); %magnitude
% BA.mag(:,2) = BVhist(~isnan(BVhist(:,4)));


%% RMSE
BA.ap_sq_diff = ( BA.ap(:,1) - BA.ap(:,2) ).^2;
BA.ap_mse     = sum( BA.ap_sq_diff ./ length(BA.ap) );
BA.ap_rmse    = sqrt( BA.ap_mse );
BA.ap_imrmse  = sqrt( immse( BA.ap(:,1), BA.ap(:,2) ) );

BA.fh_sq_diff = ( BA.fh(:,1) - BA.fh(:,2) ).^2;
BA.fh_mse     = sum( BA.fh_sq_diff ./ length(BA.fh) );
BA.fh_rmse    = sqrt( BA.fh_mse );
BA.fh_imrmse  = sqrt( immse( BA.fh(:,1), BA.fh(:,2) ) );

BA.rl_sq_diff = ( BA.rl(:,1) - BA.rl(:,2) ).^2;
BA.rl_mse     = sum( BA.rl_sq_diff ./ length(BA.rl) );
BA.rl_rmse    = sqrt( BA.rl_mse );
BA.rl_imrmse  = sqrt( immse( BA.rl(:,1), BA.rl(:,2) ) );

BA.mag_sq_diff = ( BA.mag(:,1) - BA.mag(:,2) ).^2;
BA.mag_mse     = sum( BA.mag_sq_diff ./ length(BA.mag) );
BA.mag_rmse    = sqrt( BA.mag_mse );
BA.mag_imrmse  = sqrt( immse( BA.mag(:,1), BA.mag(:,2) ) );

figure; hold on;
subplot(2,2,1);
plot(BA.ap(1:50:end,1),BA.ap(1:50:end,2),'o'); axis([-10 10 -10 10]);
subplot(2,2,2);
plot(BA.fh(1:50:end,1),BA.fh(1:50:end,2),'o'); axis([-10 10 -10 10]);
subplot(2,2,3);
plot(BA.rl(1:50:end,1),BA.rl(1:50:end,2),'o'); axis([-10 10 -10 10]);
subplot(2,2,4);
plot(BA.mag(1:50:end,1),BA.mag(1:50:end,2),'o'); axis([-10 10 -10 10]);




