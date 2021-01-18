function T_all = fcmr_4dflow_vessel_analysis( studyDir, fcmrNum, usePolyCorr, roiDir )

%% 4D Cardiac Flow - fcmr vessel analysis
%
% - 2D ROIs drawn around vessels / regions in MITK
% - requires .nii.gz files saved in fcmrDir/roi
% - ideally ROI files will be named according to vessel, ie: DAo.nii.gz
% - script uses file names for titles in figures.
%
% See also:
%   Wrapper script: ...\CODE_Projects\fetal_cmr_4dflow\wip\fcmr_4dflow_vesselanalysis_wrapper.m

% Tom Roberts (t.roberts@kcl.ac.uk)


close all
% clear

%% Admin

% fcmrNum = 191;
try
    fcmrDir = fullfile( studyDir, ['fcmr' num2str(fcmrNum)] );
    cd(fcmrDir);
catch
    fcmrDir = fullfile( studyDir, ['c_fcmr' num2str(fcmrNum)] );
    cd(fcmrDir);
end

velDir = 'vel_vol';

if nargin < 3
    usePolyCorr = false;
    roiDir = false;
elseif nargin < 4
    roiDir = false;
end


%% Load cine volume

cd('cine_vol')
if ~isfile('cine_vol-RESLICE.nii.gz')
    reslice_nii('cine_vol.nii.gz', 'cine_vol-RESLICE.nii.gz');
end
cine_nii = load_nii('cine_vol-RESLICE.nii.gz');

nFrame = size(cine_nii.img,4);


%% Load velocity volume

cd(['../' velDir]);

if usePolyCorr == false

    %%% straight out of SVRTK
    if ~isfile('velocity-final-RESLICE-0.nii.gz') || ~isfile('velocity-final-RESLICE-1.nii.gz') || ~isfile('velocity-final-RESLICE-2.nii.gz') 
            reslice_nii('velocity-final-0.nii.gz', 'velocity-final-RESLICE-0.nii.gz');
        reslice_nii('velocity-final-1.nii.gz', 'velocity-final-RESLICE-1.nii.gz');
        reslice_nii('velocity-final-2.nii.gz', 'velocity-final-RESLICE-2.nii.gz');
    end
    
    velx_nii = load_nii('velocity-final-RESLICE-0.nii.gz');
    vely_nii = load_nii('velocity-final-RESLICE-1.nii.gz');
    velz_nii = load_nii('velocity-final-RESLICE-2.nii.gz');

elseif usePolyCorr == true

    %%% velocity volume drift correction
    if ~isfile('velocity-final-polyCorr-RESLICE-0.nii.gz') || ~isfile('velocity-final-polyCorr-RESLICE-1.nii.gz') || ~isfile('velocity-final-polyCorr-RESLICE-2.nii.gz') 
        reslice_nii('velocity-final-polyCorr-0.nii.gz', 'velocity-final-polyCorr-RESLICE-0.nii.gz');
        reslice_nii('velocity-final-polyCorr-1.nii.gz', 'velocity-final-polyCorr-RESLICE-1.nii.gz');
        reslice_nii('velocity-final-polyCorr-2.nii.gz', 'velocity-final-polyCorr-RESLICE-2.nii.gz');
    end

    velx_nii = load_nii('velocity-final-polyCorr-RESLICE-0.nii.gz');
    vely_nii = load_nii('velocity-final-polyCorr-RESLICE-1.nii.gz');
    velz_nii = load_nii('velocity-final-polyCorr-RESLICE-2.nii.gz');
    
end

Vx = velx_nii.img;
Vy = vely_nii.img;
Vz = velz_nii.img;

% cm/s
Vx = 1e2 .* Vx;
Vy = 1e2 .* Vy; 
Vz = 1e2 .* Vz;


%% Load blood pool mask

cd('../mask');
maskFileName = 'mask_blood_pool';

% open mask / reslice if necessary
if ~isfile([maskFileName '-RESLICE.nii.gz'])
    reslice_nii([maskFileName '.nii.gz'], [maskFileName '-RESLICE.nii.gz']);
end

mask = load_nii([maskFileName '-RESLICE.nii.gz']);
mask.img = double(mask.img);
for ii = 1:nFrame; mask.img(:,:,:,ii) = mask.img(:,:,:,1); end


%% Load vessel ROIs

cd(fcmrDir);

% Folder check
if roiDir == false
    if exist('roi_v2') == 7
        roiDir = 'roi_v2';
        cd(roiDir);
    else
        disp('Version 2 folder does not exist ... ')
        roiDir = 'roi';
        cd(roiDir);
    end
else
    disp(['Using ROIs in folder: ' roiDir]);
    cd(roiDir);
end

dirNames = dir('*.nii.gz');
for ii = 1:numel(dirNames)
    roiNames(ii) = cellstr(dirNames(ii).name(1:end-7));
end

for ii = 1:numel(roiNames)
    ROI.(roiNames{ii}) = load_nii([roiNames{ii} '.nii.gz']);
    ROI.(roiNames{ii}).img = double(ROI.(roiNames{ii}).img);
    
    % find frame containing mask
    roiIdx    = find(ROI.(roiNames{ii}).img);
    [~,~,~,roiFrameNum] = ind2sub( size(ROI.(roiNames{ii}).img) , roiIdx );
    roiFrameNum = unique(roiFrameNum);
    if numel(roiFrameNum) > 1
        error('The mask contains ROIs in more than one frame!');
    end
    
    % propagate mask to all frames
    for nn = 1:nFrame
        ROI.(roiNames{ii}).img(:,:,:,nn) = ROI.(roiNames{ii}).img(:,:,:,roiFrameNum); % originally roiFrameNum = 1
    end 
end


%% Resize
% FIXME: for some reason mask is 1 voxel larger than original volume...? 

for tt = 1:nFrame
    Vx_re(:,:,:,tt) = imresize3(Vx(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    Vy_re(:,:,:,tt) = imresize3(Vy(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    Vz_re(:,:,:,tt) = imresize3(Vz(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    mask_re.img(:,:,:,tt) = imresize3(mask.img(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));    
    
    for ii = 1:numel(roiNames)
        ROI.([roiNames{ii} '_re']).img(:,:,:,tt) = imresize3(ROI.(roiNames{ii}).img(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    end
    
end

Vx = Vx_re; clear Vx_re;
Vy = Vy_re; clear Vy_re;
Vz = Vz_re; clear Vz_re;
mask.img = mask_re.img; clear mask_re;

for ii = 1:numel(roiNames)
    ROI.(roiNames{ii}).img = ROI.([roiNames{ii} '_re']).img; ROI = rmfield(ROI,[roiNames{ii} '_re']);
end


%% Calculate magnitude /sum velocity volume

Vsum = (Vx + Vy + Vz);
Vmag = sqrt(Vx.^2 + Vy.^2 + Vz.^2);


%% Measure velocity in ROIs

% Fetch velocities for each ROI
for ii = 1:numel(roiNames)
    ROI.(roiNames{ii}) = roi2vel( ROI.(roiNames{ii}), Vsum, Vmag, nFrame );
end    

% Automatically determine frame shift to align systole/diastole
for ii = 1:numel(roiNames)
    velMean = ROI.(roiNames{ii}).vmag;
    systoleCardiacPhaseFraction = 0.3;
    targetSystoleFrameNo        = round( systoleCardiacPhaseFraction * nFrame );
    maxVelocityFrameNo          = find( abs(velMean) == max(abs(velMean)), 1 );
    nFrameShift(ii)             = targetSystoleFrameNo - maxVelocityFrameNo;
end
   
% Keep shift consistent between all ROIs
medianFrameShift = round(median(nFrameShift));


%% Save to analysis folder 

cd( fullfile (fcmrDir , roiDir) )

if usePolyCorr == true
    mkdir('vessel_analysis_polyCorr');
    cd('vessel_analysis_polyCorr');
elseif usePolyCorr == false
    mkdir('vessel_analysis');
    cd('vessel_analysis');
end


%% Temporal information

rr = cine_nii.hdr.dime.pixdim(5) * 1e3 * nFrame; % ms
dt = cine_nii.hdr.dime.pixdim(5) * 1e3;          % ms
time = linspace(0,rr,nFrame);                    % ms

fileID = fopen(['fcmr' num2str(fcmrNum) '_frame_times.txt'],'w');
fprintf(fileID,'Median frame shift: \t');
fprintf(fileID,'%.4f', medianFrameShift);
fprintf(fileID,'\n\n');
fprintf(fileID,'Frame times (ms): \n');
fprintf(fileID,'%.4f\t', time);
fclose(fileID);


% Make subplots
%% \\ VELOCITY SUM

hFig = figure;
ftSz = 14;

fileID = fopen(['fcmr' num2str(fcmrNum) '_velocity_sum_vessels.txt'],'w');
fprintf(fileID,[ 'fcmr' num2str(fcmrNum) ': \n']);
fprintf(fileID,'Mean / STD \n\n');

for ii = 1:numel(roiNames)

    subplot(2,4,ii);
    velMean = circshift(ROI.(roiNames{ii}).vsum,medianFrameShift);
    velStd  = circshift(ROI.(roiNames{ii}).vstd,medianFrameShift);
    velMax  = circshift(ROI.(roiNames{ii}).vmax,medianFrameShift);
    velMin  = circshift(ROI.(roiNames{ii}).vmin,medianFrameShift);
    
%     % flip graphs if 80% of values negative
%     if numel(find(velMean < 0)) > nFrame * 0.8
%         velMean = velMean .* -1;
%         velStd  = velStd .* -1;
%     end
    
    c = [ 0.7 0.8 0.9 ] ;
    x = [1:nFrame,flip(1:nFrame,2)];
    y = [velMean+velStd,flip(velMean-velStd,2)]; % STDEV bounds
%     y = [velMax,flip(velMin,2)]; % max/min bounds
    patch(x,y,c.^0.5,'EdgeColor',c,'LineWidth',0.1),
    line( time, zeros(size(1:numel(time))), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5 ), 
    line( 1:nFrame, velMean, 'Color', 'b', 'LineWidth', 2 ), 
    xlabel( 'Frame number' , 'FontSize', ftSz ), 
    ylabel( 'Velocity (cm/s)' , 'FontSize', ftSz ),
    title( roiNames{ii} , 'FontSize', ftSz),
    axis('square');
    hAx = gca; hAx.XLim = [1 25]; hAx.YLim = [-5 30]; hAx.FontSize = ftSz; hAx.Box = 'on';
    
    fprintf(fileID,[ roiNames{ii} ': \n']);
    fprintf(fileID,'%.4f\t', velMean);
    fprintf(fileID,'%.4f\t', velStd);
    fprintf(fileID,'\n\n');
    
    % collate for spreadsheet
    VSum.(roiNames{ii}).velMean = velMean;
    VSum.(roiNames{ii}).velStd = velStd;
    
end

screenSize = get( 0, 'screensize' );
hFig.Position(1) = 0.05 * screenSize(3);
hFig.Position(3) = 0.9 * screenSize(3);
hFig.Position(2) = 0.1 * screenSize(4);
hFig.Position(4) = 0.8 * screenSize(4);

saveas(gcf,['fcmr' num2str(fcmrNum) '_velocity_sum_vessels.fig']);
saveas(gcf,['fcmr' num2str(fcmrNum) '_velocity_sum_vessels.png']);

fclose(fileID);


%% \\ VELOCITY MAGNITUDE

hFig = figure;
ftSz = 14;

fileID = fopen(['fcmr' num2str(fcmrNum) '_velocity_magnitude_vessels.txt'],'w');
fprintf(fileID,[ 'fcmr' num2str(fcmrNum) ': \n']);
fprintf(fileID,'Mean / STD \n\n');
    
for ii = 1:numel(roiNames)

    subplot(2,4,ii);
    velMean = circshift(ROI.(roiNames{ii}).vmag,medianFrameShift);
    velStd  = circshift(ROI.(roiNames{ii}).vstd,medianFrameShift);
    velMax  = circshift(ROI.(roiNames{ii}).vmax,medianFrameShift);
    velMin  = circshift(ROI.(roiNames{ii}).vmin,medianFrameShift);
    
    c = [ 0.7 0.8 0.9 ] ;
    x = [1:nFrame,flip(1:nFrame,2)];
    y = [velMean+velStd,flip(velMean-velStd,2)]; % STDEV bounds
%     y = [velMax,flip(velMin,2)]; % max/min bounds
    patch(x,y,c.^0.5,'EdgeColor',c,'LineWidth',0.1),
    line( time, zeros(size(1:numel(time))), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5 ), 
    line( 1:nFrame, velMean, 'Color', 'b', 'LineWidth', 2 ), 
    xlabel( 'Frame number' , 'FontSize', ftSz ), 
    ylabel( 'Velocity Magnitude (cm/s)' , 'FontSize', ftSz ),
    title( roiNames{ii} , 'FontSize', ftSz),
    axis('square');
    hAx = gca; hAx.XLim = [1 25]; hAx.YLim = [-5 30]; hAx.FontSize = ftSz; hAx.Box = 'on';
    
    fprintf(fileID,[ roiNames{ii} ': \n']);
    fprintf(fileID,'%.4f\t', velMean);
    fprintf(fileID,'%.4f\t', velStd);
    fprintf(fileID,'\n\n');
    
    % collate for spreadsheet
    VMag.(roiNames{ii}).velMean = velMean;
    VMag.(roiNames{ii}).velStd = velStd;
    
end

screenSize = get( 0, 'screensize' );
hFig.Position(1) = 0.05 * screenSize(3);
hFig.Position(3) = 0.9 * screenSize(3);
hFig.Position(2) = 0.1 * screenSize(4);
hFig.Position(4) = 0.8 * screenSize(4);

saveas(gcf,['fcmr' num2str(fcmrNum) '_velocity_magnitude_vessels.fig']);
saveas(gcf,['fcmr' num2str(fcmrNum) '_velocity_magnitude_vessels.png']);

% close;

fclose(fileID);


%% \\ FLOW MAGNITUDE

% get voxel size
voxRes = cine_nii.hdr.dime.pixdim(2:4);
voxArea = voxRes(1) * voxRes(2); % mm^2

hFig = figure;
ftSz = 14;

fileID = fopen(['fcmr' num2str(fcmrNum) '_flow_magnitude_vessels.txt'],'w');
fprintf(fileID,[ 'fcmr' num2str(fcmrNum) ': \n']);
fprintf(fileID,'Mean / STD / Total Mean / Total STD / Stroke Volume \n\n');

for ii = 1:numel(roiNames)

    velMean = circshift(ROI.(roiNames{ii}).vmag,medianFrameShift);
    velStd  = circshift(ROI.(roiNames{ii}).vstd,medianFrameShift);
    velMax  = circshift(ROI.(roiNames{ii}).vmax,medianFrameShift);
    velMin  = circshift(ROI.(roiNames{ii}).vmin,medianFrameShift);
    
    % convert current velocity through ROI to flow
    roiMask = ROI.(roiNames{ii}).img(:,:,:,1); % mask identical between frames
    roiArea = voxArea * nnz( roiMask ) * ones( size( velMean ) );  % mm^2
    vel2flow = @( v, a ) ( v ) .* ( a * 1e-2 );  % ml/s
    flowMean = vel2flow( velMean, roiArea );
    flowStd = vel2flow( velStd, roiArea );
    flowMax = vel2flow( velMax, roiArea );
    flowMin = vel2flow( velMin, roiArea );
    flowTotalMean = mean( flowMean ) * 60; % ml/min
    flowTotalStd = mean( flowStd ) * 60; % ml/min   
    strokeVolume = sum( flowMean * 1e-3 .* dt );  % ml   
    
    subplot(2,4,ii);
    c = [ 0.7 0.8 0.9 ] ;

    x = [time,flip(time,2)];
    y = [flowMean+flowStd,flip(flowMean-flowStd,2)]; % STDEV bounds
%     y = [flowMax,flip(flowMin,2)]; % max/min bounds
    patch(x,y,c.^0.5,'EdgeColor',c,'LineWidth',0.1),
    line( time, zeros(size(1:numel(time))), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5 ), 
    line( time, flowMean, 'Color', 'b', 'LineWidth', 2 ), 
    xlabel( 'Time (ms)' , 'FontSize', ftSz ), 
    
%     % frame number instead:
%     x = [1:nFrame,flip(1:nFrame,2)];
%     line( 1:nFrame, flowMean, 'Color', 'b', 'LineWidth', 2 ), 
%     xlabel( 'Frame number' , 'FontSize', ftSz ), 
    
    ylabel( 'Flow Magnitude (ml/s)' , 'FontSize', ftSz ),
    title( roiNames{ii} , 'FontSize', ftSz),
    axis('square');
    hAx = gca; hAx.XLim = [0 rr]; hAx.YLim = [-5 15]; hAx.FontSize = ftSz; hAx.Box = 'on';
    
    txt = {['Mean Flow = ' num2str(round(flowTotalMean)) ' ml/min' ]}; %['Stroke Volume = ' num2str(round(strokeVolume)) ' ml' ] 
    text(0.05.*hAx.XLim(2),0.95.*hAx.YLim(2),txt);
    
    disp(['Mean Flow in ' roiNames{ii} ' = ' num2str(round(flowTotalMean)) ' ml/min' ]);
    
    fprintf(fileID,[ roiNames{ii} ': \n']);
    fprintf(fileID,'%.4f\t', flowMean);
    fprintf(fileID,'%.4f\t', flowStd);
    fprintf(fileID,'%.4f\n', flowTotalMean);
    fprintf(fileID,'%.4f\n', flowTotalStd);
    fprintf(fileID,'%.4f\n', strokeVolume);
    fprintf(fileID,'\n\n');
    
    % collate for spreadsheet
    FMag.(roiNames{ii}).flowMean = flowMean;
    FMag.(roiNames{ii}).flowStd = flowStd;
    FMag.(roiNames{ii}).flowTotalMean = flowTotalMean;
    
end

screenSize = get( 0, 'screensize' );
hFig.Position(1) = 0.05 * screenSize(3);
hFig.Position(3) = 0.9 * screenSize(3);
hFig.Position(2) = 0.1 * screenSize(4);
hFig.Position(4) = 0.8 * screenSize(4);

saveas(gcf,['fcmr' num2str(fcmrNum) '_flow_magnitude_vessels.fig']);
saveas(gcf,['fcmr' num2str(fcmrNum) '_flow_magnitude_vessels.png']);

fclose(fileID);


%% \\ FLOW SUM

% get voxel size
voxRes = cine_nii.hdr.dime.pixdim(2:4);
voxArea = voxRes(1) * voxRes(2); % mm^2

hFig = figure;
ftSz = 14;

fileID = fopen(['fcmr' num2str(fcmrNum) '_flow_sum_vessels.txt'],'w');
fprintf(fileID,[ 'fcmr' num2str(fcmrNum) ': \n']);
fprintf(fileID,'Mean / STD / Total Mean / Total STD / Stroke Volume \n\n');

for ii = 1:numel(roiNames)

    velMean = circshift(ROI.(roiNames{ii}).vsum,medianFrameShift);
    velStd  = circshift(ROI.(roiNames{ii}).vstd,medianFrameShift);
    velMax  = circshift(ROI.(roiNames{ii}).vmax,medianFrameShift);
    velMin  = circshift(ROI.(roiNames{ii}).vmin,medianFrameShift);
    
    % convert current velocity through ROI to flow
    roiMask = ROI.(roiNames{ii}).img(:,:,:,1); % mask identical between frames
    roiArea = voxArea * nnz( roiMask ) * ones( size( velMean ) );  % mm^2
    vel2flow = @( v, a ) ( v ) .* ( a * 1e-2 );  % ml/s
    flowMean = vel2flow( velMean, roiArea );
    flowStd = vel2flow( velStd, roiArea );
    flowMax = vel2flow( velMax, roiArea );
    flowMin = vel2flow( velMin, roiArea );
    flowTotalMean = mean( flowMean ) * 60; % ml/min
    flowTotalStd = mean( flowStd ) * 60; % ml/min   
    strokeVolume = sum( flowMean * 1e-3 .* dt );  % ml   
    
    subplot(2,4,ii);
    c = [ 0.7 0.8 0.9 ] ;

    x = [time,flip(time,2)];
    y = [flowMean+flowStd,flip(flowMean-flowStd,2)]; % STDEV bounds
%     y = [flowMax,flip(flowMin,2)]; % max/min bounds
    patch(x,y,c.^0.5,'EdgeColor',c,'LineWidth',0.1),
    line( time, zeros(size(1:numel(time))), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5 ), 
    line( time, flowMean, 'Color', 'b', 'LineWidth', 2 ), 
    xlabel( 'Time (ms)' , 'FontSize', ftSz ), 
    
%     % frame number instead:
%     x = [1:nFrame,flip(1:nFrame,2)];
%     line( 1:nFrame, flowMean, 'Color', 'b', 'LineWidth', 2 ), 
%     xlabel( 'Frame number' , 'FontSize', ftSz ), 
    
    ylabel( 'Flow (ml/s)' , 'FontSize', ftSz ),
    title( roiNames{ii} , 'FontSize', ftSz),
    axis('square');
    hAx = gca; hAx.XLim = [0 rr]; hAx.YLim = [-5 15]; hAx.FontSize = ftSz; hAx.Box = 'on';
    
    txt = {['Mean Flow = ' num2str(round(flowTotalMean)) ' ml/min' ]}; %['Stroke Volume = ' num2str(round(strokeVolume)) ' ml' ] 
    text(0.05.*hAx.XLim(2),0.95.*hAx.YLim(2),txt);
    
    disp(['Mean Flow in ' roiNames{ii} ' = ' num2str(round(flowTotalMean)) ' ml/min' ]);
    
    fprintf(fileID,[ roiNames{ii} ': \n']);
    fprintf(fileID,'%.4f\t', flowMean);
    fprintf(fileID,'%.4f\t', flowStd);
    fprintf(fileID,'%.4f\n', flowTotalMean);
    fprintf(fileID,'%.4f\n', flowTotalStd);
    fprintf(fileID,'%.4f\n', strokeVolume);
    fprintf(fileID,'\n\n');
    
    % collate for spreadsheet
    FSum.(roiNames{ii}).flowMean = flowMean;
    FSum.(roiNames{ii}).flowStd = flowStd;
    FSum.(roiNames{ii}).flowTotalMean = flowTotalMean;
    
end

screenSize = get( 0, 'screensize' );
hFig.Position(1) = 0.05 * screenSize(3);
hFig.Position(3) = 0.9 * screenSize(3);
hFig.Position(2) = 0.1 * screenSize(4);
hFig.Position(4) = 0.8 * screenSize(4);

saveas(gcf,['fcmr' num2str(fcmrNum) '_flow_sum_vessels.fig']);
saveas(gcf,['fcmr' num2str(fcmrNum) '_flow_sum_vessels.png']);

fclose(fileID);


%% Save spreadsheet
% cd(fcmrDir);
% cd(roiDir);

ssFileName = ['fcmr' num2str(fcmrNum) '_flow_velocity_values.xlsx'];

% save ROI per sheet
for ii = 1:numel(roiNames)
    
    T = table( time', ...
               VSum.(roiNames{ii}).velMean', ...
               VSum.(roiNames{ii}).velStd', ...
               VMag.(roiNames{ii}).velMean', ...
               VMag.(roiNames{ii}).velStd', ...
               FSum.(roiNames{ii}).flowMean', ...
               FSum.(roiNames{ii}).flowStd', ...
               FSum.(roiNames{ii}).flowTotalMean .* ones(nFrame,1), ...
               FMag.(roiNames{ii}).flowMean', ...
               FMag.(roiNames{ii}).flowStd', ...
               FMag.(roiNames{ii}).flowTotalMean .* ones(nFrame,1) ...
               );
           
    T.Properties.VariableNames = {'FrameTimes' , ...
                                  ['VSum_' roiNames{ii} '_Mean'], ...
                                  ['VSum_' roiNames{ii} '_STD'], ...
                                  ['VMag_' roiNames{ii} '_Mean'], ...
                                  ['VMag_' roiNames{ii} '_STD'], ...
                                  ['FSum_' roiNames{ii} '_Mean'], ...
                                  ['FSum_' roiNames{ii} '_STD'], ...
                                  ['FSum_' roiNames{ii} '_MeanTotal'], ...
                                  ['FMag_' roiNames{ii} '_Mean'], ...
                                  ['FMag_' roiNames{ii} '_STD'], ...
                                  ['FMag_' roiNames{ii} '_MeanTotal'] ...
                                  };
    
    writetable(T,ssFileName,'Sheet',roiNames{ii});
    warning('off','last'); % turn off sheet created warning    
end


% save vessels per sheet
% REALLY UGLY PROGRAMMING... argh.

analysisTypes = {'VSum','VMag','FSum','FMag'};

for aa = 1:numel(analysisTypes)
    
    switch analysisTypes{aa}
        case 'VSum'
            for ii = 1:numel(roiNames)
                vesselsMatrix(:,ii) = VSum.(roiNames{ii}).velMean';
                vesselsMatrix(:,numel(roiNames)+ii) = VSum.(roiNames{ii}).velStd';
            end
        case 'VMag'
            for ii = 1:numel(roiNames)
                vesselsMatrix(:,ii) = VMag.(roiNames{ii}).velMean';
                vesselsMatrix(:,numel(roiNames)+ii) = VMag.(roiNames{ii}).velStd';
            end
        case 'FSum'
            for ii = 1:numel(roiNames)
                vesselsMatrix(:,ii) = FSum.(roiNames{ii}).flowMean';
                vesselsMatrix(:,numel(roiNames)+ii) = FSum.(roiNames{ii}).flowStd';
            end
        case 'FMag'
            for ii = 1:numel(roiNames)
                vesselsMatrix(:,ii) = FMag.(roiNames{ii}).flowMean';
                vesselsMatrix(:,numel(roiNames)+ii) = FMag.(roiNames{ii}).flowStd';
            end
    end

    T = array2table([time', vesselsMatrix]);
    T.Properties.VariableNames{1} = 'FrameTimes';

    for ii = 1:numel(roiNames)
        T.Properties.VariableNames{ii+1} = [roiNames{ii} '_Mean'];
        T.Properties.VariableNames{numel(roiNames)+ii+1} = [roiNames{ii} '_STD'];
    end

    writetable(T,ssFileName,'Sheet',analysisTypes{aa});
    warning('off','last'); % turn off sheet created warning
    
    % save tables for function output
    switch analysisTypes{aa}
        case 'VSum'
            T_all.VSUM = T;
        case 'VMag'
            T_all.VMAG = T;
        case 'FSum'
            T_all.FSUM = T;
        case 'FMag'
            T_all.FMAG = T;
    end
end
                              

% delete Sheet1/Sheet2/Sheet3
sheetName = 'Sheet';
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(fullfile(pwd, ssFileName));
try
      % Throws an error if the sheets do not exist.
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
catch
      % Do nothing.
end
% Save, close and clean up.
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;

% fn end
end


%% subfn
function ROI = roi2vel(ROI, Vsum, Vmag, nFrame)

for tt = 1:nFrame
    ROI.vsum(tt) = mean(nonzeros(Vsum(:,:,:,tt) .* ROI.img(:,:,:,tt)));
    ROI.vmag(tt) = mean(nonzeros(Vmag(:,:,:,tt) .* ROI.img(:,:,:,tt)));
    ROI.vstd(tt) = std(nonzeros(Vmag(:,:,:,tt) .* ROI.img(:,:,:,tt)));
    ROI.vmax(tt) = max(nonzeros(Vmag(:,:,:,tt) .* ROI.img(:,:,:,tt)));
    ROI.vmin(tt) = min(nonzeros(Vmag(:,:,:,tt) .* ROI.img(:,:,:,tt)));
end

% rectify vsum
ROI.vsum = sign( sum( ROI.vsum(:) ) ) * ROI.vsum;


end %fcmr_4dflow_vessel_analysis(...)
