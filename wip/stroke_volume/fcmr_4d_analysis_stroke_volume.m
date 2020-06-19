
%% fcmr_4D_stroke_volume_analysis
%
% Created for Chloe in April 2020
%
% Aim: scripts for analysing stroke volume in Chloe's data
%   - LV volume measured at EDV and ESV in cine_vol.nii.gz
%   - Aortic flow measured in AAo using velocity cine volume
%
%

%% Setup
% fcmrIDs  = [191 194 197 201 202 213 214 242 252 256 266 305];
fcmrIDs  = 256;

for ff = fcmrIDs
    
    fcmrNum = ff;
    disp(['Running fcmr ' num2str(fcmrNum) ' ...']);
    
    % File location check
    if fcmrNum == 191 || fcmrNum == 194 || fcmrNum == 197 || fcmrNum == 201 || fcmrNum == 202 || fcmrNum == 214   
        studyDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper';
        reconDir = ['E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr' num2str(fcmrNum)];   
    elseif fcmrNum == 213
        studyDir = 'I:\fcmr_4d_chloe_recons';
        reconDir = ['I:\fcmr_4d_chloe_recons\fcmr' num2str(fcmrNum)];
    else
        studyDir = 'I:\fcmr_4d_chloe_recons';
        reconDir = ['I:\fcmr_4d_chloe_recons\c_fcmr' num2str(fcmrNum)];
    end
    cd(reconDir);
    
    % Folder/File names
    cineVolDir   = 'cine_vol';
    velVolDir    = 'vel_vol';
    strokeVolDir = 'stroke_vol';
    roilvDir     = 'roi_LV';
    
    cineVolFilename = 'cine_vol-RESLICE.nii.gz';
    lvedvFilename   = 'LVEDV.nii.gz';
    lvesvFilename   = 'LVESV.nii.gz';
    aaoedvFilename  = 'AAo_EDV.nii.gz';
    aaoesvFilename  = 'AAo_ESV.nii.gz';
    
    % use polynomial corrected data:
    usePolyCorr = true;
    
    
    %% Load cine_vol
    % cd( [reconDir '\' cineVolDir] ); % Another way of navigating directories
    cd( fullfile( reconDir, cineVolDir ) );
    nii     = load_untouch_nii( cineVolFilename );
    pixdim  = nii.hdr.dime.pixdim(2:4);          % mm
    pixvol  = pixdim(1) * pixdim(2) * pixdim(3); % mm^3
    cineVol = double(nii.img);
    
    
    %% Load LV segmentations
    cd( fullfile( reconDir, strokeVolDir, roilvDir ) );
    nii      = load_untouch_nii( lvedvFilename );
    lvedvVol = double(nii.img);
    nii      = load_untouch_nii( lvesvFilename );
    lvesvVol = double(nii.img);
    
    
    %% Load AAo segmentations
    cd( fullfile( reconDir, strokeVolDir ) );
    nii       = load_untouch_nii( aaoedvFilename );
    aaoedvVol = double(nii.img);
%     nii       = load_untouch_nii( aaoesvFilename );
%     aaoesvVol = double(nii.img);
    
    
    %% Stroke volume - volumetric calculation
    % 1 mm^3 = 0.001 mL
    volEDV = cineVol .* lvedvVol;
    numPixelsEDV = numel( find( volEDV ) );
    LVEDV = numPixelsEDV .* pixvol * 1e-3; % mL
    
    volESV = cineVol .* lvesvVol;
    numPixelsESV = numel( find( volESV ) );
    LVESV = numPixelsESV .* pixvol * 1e-3; % mL
    
    LVSV = LVEDV - LVESV;                  % LV stroke volume - mL
    LVEF = 100 * (LVSV / LVEDV);           % LV ejection fraction - %
    
    
    %% Stroke volume - flow calculation
    T_all = fcmr_4dflow_vessel_analysis( studyDir, fcmrNum, usePolyCorr, strokeVolDir );
    close all;
    
    frameTimes = T_all.FMAG{:,1};
    flowCurve  = T_all.FMAG{:,2};
    
    dt = max( frameTimes ) / numel( frameTimes ) * 1e-3; % seconds
    
    cumVol = cumtrapz( flowCurve ) * dt;
    flowSV = cumVol(end); % mL
    
    
    %% Plot flow curve
    figure;
    plot( frameTimes, flowCurve / 60, '-bo' ); % / 60 to go from mL/min to mL/second
    xlabel('Time [ms]');
    ylabel('Flow [mL/second]');
    title(['Flow curve, patient number = ' num2str(fcmrNum) ]);
    annotation('textbox', [0.4, 0.8, 0.1, 0.1], ...
        'String', ['LVSV = ' num2str(round(LVSV,2)) ' mL // 4D Flow SV = ' num2str(round(flowSV,2)) ' mL' ] );
    
    
    %% Save Stroke Volume Values
    cd( fullfile( reconDir, strokeVolDir ) );
    save( 'stroke_volume_analysis.mat', 'LVSV', 'LVEF', 'flowSV', 'frameTimes', 'flowCurve', 'dt' );
    
    disp(['Completed fcmr ' num2str(fcmrNum) '.']);
    % end for loop
end


disp('ALL COMPLETE.');























