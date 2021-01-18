

clear all; close all;


fcmrNum = 214;

%%%%% NB: need to resample vol_vol_orig for this fcmr194 - for some reason
%%%%% different imsize to cine_vol_orig. Used:
%%%%% mirtk resample-image vel_vol_0_orig.nii.gz vel_vol_0_orig_resampled.nii.gz -imsize 65 73 64



velVolRegDir      = ['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum) '_hrh_fullRecon\4dvol_reg'];
roiDir = [velVolRegDir '/roi'];

usePolyCorr = true;


%% Load registered magnitude volume
% These are -RESLICE / polycorr volumes
cd(velVolRegDir);
cine_nii_orig = load_nii('cine_vol_orig.nii.gz');
cine_nii_hhh_reg = load_nii('cine_vol_hhh_reg.nii.gz');
nFrame = size(cine_nii_orig.img,4);


%% Load registered velocity volumes
velx_nii_orig = load_nii('vel_vol_0_orig.nii.gz');
vely_nii_orig = load_nii('vel_vol_1_orig.nii.gz');
velz_nii_orig = load_nii('vel_vol_2_orig.nii.gz');

velx_nii_hhh_reg = load_nii('vel_vol_0_hhh_reg.nii.gz');
vely_nii_hhh_reg = load_nii('vel_vol_1_hhh_reg.nii.gz');
velz_nii_hhh_reg = load_nii('vel_vol_2_hhh_reg.nii.gz');

Vx_orig = velx_nii_orig.img;
Vy_orig = vely_nii_orig.img;
Vz_orig = velz_nii_orig.img;

Vx_hhh_reg = velx_nii_hhh_reg.img;
Vy_hhh_reg = vely_nii_hhh_reg.img;
Vz_hhh_reg = velz_nii_hhh_reg.img;

% cm/s
Vx_orig = 1e2 .* Vx_orig;
Vy_orig = 1e2 .* Vy_orig;
Vz_orig = 1e2 .* Vz_orig;

Vx_hhh_reg = 1e2 .* Vx_hhh_reg;
Vy_hhh_reg = 1e2 .* Vy_hhh_reg;
Vz_hhh_reg = 1e2 .* Vz_hhh_reg;


% Calculate magnitude /sum velocity volume
Vmag_orig = sqrt(Vx_orig.^2 + Vy_orig.^2 + Vz_orig.^2);
Vmag_hhh_reg = sqrt(Vx_hhh_reg.^2 + Vy_hhh_reg.^2 + Vz_hhh_reg.^2);



%% Load vessel ROIs

cd(roiDir);

% dirNames = dir('*_reg.nii.gz');
dirNames = dir('*.nii.gz');
for ii = 1:numel(dirNames)
    roiNames(ii) = cellstr(dirNames(ii).name(1:end-7));
end

for ii = 1:numel(roiNames)
    ROI.(roiNames{ii}) = load_nii([roiNames{ii} '.nii.gz']);
    ROI.(roiNames{ii}).img = double(ROI.(roiNames{ii}).img);
    
    % propagate mask to all frames
    for nn = 1:nFrame
        ROI.(roiNames{ii}).img(:,:,:,nn) = ROI.(roiNames{ii}).img(:,:,:,1);
    end
end

cd(velVolRegDir);


%% Calculate velocities

numROIs = numel(dirNames);

for ii = 1:numROIs
    
    for ss = 1:nFrame
        % original
        ROI.(roiNames{ii}).orig.Vmag_mean(:,ss) = mean( nonzeros( Vmag_orig(:,:,:,ss) .* ROI.(roiNames{ii}).img(:,:,:,ss) ) ,1 );
        
        % hhh_reg
        ROI.(roiNames{ii}).hhh_reg.Vmag_mean(:,ss) = mean( nonzeros( Vmag_hhh_reg(:,:,:,ss) .* ROI.(roiNames{ii}).img(:,:,:,ss) ) ,1 );
    end
    
end


%% Calculate Flow Rates

voxRes = cine_nii_orig.hdr.dime.pixdim(2:4);
voxArea = voxRes(1) * voxRes(2); % mm^2
% vel2flow = @( v, a ) ( v ) .* ( a * 1e-2 );  % ml/s

for ii = 1:numROIs
    
    % num voxels in ROI
    ROI.(roiNames{ii}).numVoxels = numel(nonzeros(ROI.(roiNames{ii}).img(:,:,:,1)));
    
    % ROI area = voxArea * numVoxels --- mm^2
    ROI.(roiNames{ii}).roiArea = voxArea * ROI.(roiNames{ii}).numVoxels;
    
end

% Flow = vel * area
for ii = 1:numROIs
    
    for ss = 1:nFrame
        % original
        ROI.(roiNames{ii}).orig.Fmag_mean(:,ss) = ROI.(roiNames{ii}).orig.Vmag_mean(:,ss) * ROI.(roiNames{ii}).roiArea * (60 * 1e-2); % ml/min
        
        % hhh_reg

        ROI.(roiNames{ii}).hhh_reg.Fmag_mean(:,ss) = ROI.(roiNames{ii}).hhh_reg.Vmag_mean(:,ss) * ROI.(roiNames{ii}).roiArea * (60 * 1e-2); % ml/min  
    end
    
end



%% Plot

% Velocity

figure;
for ii = 1:numROIs
    
    subplot(2,3,ii)
    hold on;
    plot(ROI.(roiNames{ii}).orig.Vmag_mean,'o-');
    plot(ROI.(roiNames{ii}).hhh_reg.Vmag_mean,'ro-');
    title(roiNames{ii});
    legend('Standard k-t recon','HRH k-t recon');
    xlabel('Frame no.');
    ylabel('Velocity [cm/s]');
    
end

set(gcf, 'Position', [1, 31, 1920, 1093]);

saveas(gcf,'standard_kt_v_HRH_kt_velocity_plots.fig');
saveas(gcf,'standard_kt_v_HRH_kt_velocity_plots.png');

% Flow Rate

figure;
for ii = 1:numROIs
    
    subplot(2,3,ii)
    hold on;
    plot(ROI.(roiNames{ii}).orig.Fmag_mean,'o-');
    plot(ROI.(roiNames{ii}).hhh_reg.Fmag_mean,'ro-');
    title(roiNames{ii});
    legend('Standard k-t recon','HRH k-t recon');
    xlabel('Frame no.');
    ylabel('Flow rate [ml/min]');
    
end

set(gcf, 'Position', [1, 31, 1920, 1093]);

saveas(gcf,'standard_kt_v_HRH_kt_flow_plots.fig');
saveas(gcf,'standard_kt_v_HRH_kt_flow_plots.png');



%% Mean / STD

for ii = 1:numROIs
    
    % Mean
    ROI.(roiNames{ii}).orig.Vmag_totalVelMean = mean(ROI.(roiNames{ii}).orig.Vmag_mean);
    ROI.(roiNames{ii}).hhh_reg.Vmag_totalVelMean = mean(ROI.(roiNames{ii}).hhh_reg.Vmag_mean);

    ROI.(roiNames{ii}).orig.Fmag_totalFlowMean = mean(ROI.(roiNames{ii}).orig.Fmag_mean);
    ROI.(roiNames{ii}).hhh_reg.Fmag_totalFlowMean = mean(ROI.(roiNames{ii}).hhh_reg.Fmag_mean);
    
    
    % STD
    ROI.(roiNames{ii}).orig.Vmag_totalVelStd = std(ROI.(roiNames{ii}).orig.Vmag_mean);
    ROI.(roiNames{ii}).hhh_reg.Vmag_totalVelStd = std(ROI.(roiNames{ii}).hhh_reg.Vmag_mean);
    
    ROI.(roiNames{ii}).orig.Fmag_totalFlowStd = std(ROI.(roiNames{ii}).orig.Fmag_mean);
    ROI.(roiNames{ii}).hhh_reg.Fmag_totalFlowStd = std(ROI.(roiNames{ii}).hhh_reg.Fmag_mean);
    
end



for ii = 1:numROIs
    Vessel{ii,1} = roiNames{ii};
    
    ORIG_Flow_Mean(ii,1) = ROI.(roiNames{ii}).orig.Fmag_totalFlowMean;
    ORIG_Flow_Stdev(ii,1) = ROI.(roiNames{ii}).orig.Fmag_totalFlowStd;
    HRH_Flow_Mean(ii,1) = ROI.(roiNames{ii}).hhh_reg.Fmag_totalFlowMean;
    HRH_Flow_Stdev(ii,1) = ROI.(roiNames{ii}).hhh_reg.Fmag_totalFlowStd;
    
    ORIG_Velocity_Mean(ii,1) = ROI.(roiNames{ii}).orig.Vmag_totalVelMean;
    ORIG_Velocity_Stdev(ii,1) = ROI.(roiNames{ii}).orig.Vmag_totalVelStd;
    HRH_Velocity_Mean(ii,1) = ROI.(roiNames{ii}).hhh_reg.Vmag_totalVelMean;
    HRH_Velocity_Stdev(ii,1) = ROI.(roiNames{ii}).hhh_reg.Vmag_totalVelStd;
end

velResultsTable = table( cell2table(Vessel), ...
    ORIG_Flow_Mean, HRH_Flow_Mean, ORIG_Flow_Stdev, HRH_Flow_Stdev, ...
    ORIG_Velocity_Mean, HRH_Velocity_Mean, ORIG_Velocity_Stdev, HRH_Velocity_Stdev );


%% Save

save('velocity_analysis.mat','velResultsTable');





