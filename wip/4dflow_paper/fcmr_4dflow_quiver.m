%% fcmr_4D_velocity_quiver
%
%
%


%% Admin
fNum = 194;
fcmrDir = ['C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr' num2str(fNum)];
cd(fcmrDir);

velDir = 'vel_vol';


%% Load cine volume
cd('cine_vol')
if ~isfile('cine_vol-RESLICE.nii.gz')
    reslice_nii('cine_vol.nii.gz', 'cine_vol-RESLICE.nii.gz');
end
cine_nii = load_nii('cine_vol-RESLICE.nii.gz');

nFrame = size(cine_nii.img,4);


%% Load velocity volume
cd(['../' velDir]);

% %%% straight out of SVRTK
% reslice_nii('velocity-final-0.nii.gz', 'velocity-final-RESLICE-0.nii.gz');
% reslice_nii('velocity-final-1.nii.gz', 'velocity-final-RESLICE-1.nii.gz');
% reslice_nii('velocity-final-2.nii.gz', 'velocity-final-RESLICE-2.nii.gz');
% 
% velx_nii = load_nii('velocity-final-RESLICE-0.nii.gz');
% vely_nii = load_nii('velocity-final-RESLICE-1.nii.gz');
% velz_nii = load_nii('velocity-final-RESLICE-2.nii.gz');


%%% velocity volume polyCorr version
if ~isfile('velocity-final-polyCorr-RESLICE-0.nii.gz') || ~isfile('velocity-final-polyCorr-RESLICE-1.nii.gz') || ~isfile('velocity-final-polyCorr-RESLICE-2.nii.gz') 
    reslice_nii('velocity-final-polyCorr-0.nii.gz', 'velocity-final-polyCorr-RESLICE-0.nii.gz');
    reslice_nii('velocity-final-polyCorr-1.nii.gz', 'velocity-final-polyCorr-RESLICE-1.nii.gz');
    reslice_nii('velocity-final-polyCorr-2.nii.gz', 'velocity-final-polyCorr-RESLICE-2.nii.gz');
end

velx_nii = load_nii('velocity-final-polyCorr-RESLICE-0.nii.gz');
vely_nii = load_nii('velocity-final-polyCorr-RESLICE-1.nii.gz');
velz_nii = load_nii('velocity-final-polyCorr-RESLICE-2.nii.gz');


%% Load blood pool mask

% BP mask automatically created in fcmr_pc_velocity_correction.m

cd('../mask');

if ~isfile('mask_blood_pool-RESLICE.nii.gz')
    reslice_nii('mask_blood_pool.nii.gz', 'mask_blood_pool-RESLICE.nii.gz');
end

mask_blood_pool = load_nii('mask_blood_pool-RESLICE.nii.gz');
mask_blood_pool.img = double(mask_blood_pool.img);
% se = strel('cube',7); mask_sphere_eroded_nii.img = imerode(mask_blood_pool.img, se);


%% Make component velocity volumes for VTK
Vx = velx_nii.img;
Vy = vely_nii.img;
Vz = velz_nii.img;

% cm/s
Vx = 1e2 .* Vx;
Vy = 1e2 .* Vy; 
Vz = 1e2 .* Vz;


%% Resize
% FIXME: for some reason mask is 1 voxel larger than original volume...? 
for tt = 1:nFrame
    Vx_re(:,:,:,tt) = imresize3(Vx(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    Vy_re(:,:,:,tt) = imresize3(Vy(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    Vz_re(:,:,:,tt) = imresize3(Vz(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    mask_blood_pool_re.img(:,:,:,tt) = imresize3(mask_blood_pool.img(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
end

Vx = Vx_re; clear Vx_re;
Vy = Vy_re; clear Vy_re;
Vz = Vz_re; clear Vz_re;
mask_blood_pool.img = mask_blood_pool_re.img; clear mask_blood_pool_re;


%% Crop phantom volume using blood pool mask

% magnitude
cine_masked_nii.img = double(cine_nii.img) .* mask_blood_pool.img;

% velocity
Vx_masked = Vx .* mask_blood_pool.img;
Vy_masked = Vy .* mask_blood_pool.img;
Vz_masked = Vz .* mask_blood_pool.img;

Vmag_masked = sqrt(Vx_masked.^2 + Vy_masked.^2 + Vz_masked.^2);

% cm/s
% nb: QFlow already in cm/s


%% Multiply by -1 to match MRtrix
%%%%%%%%%%%%%%%%%%%%%%%
%%% SUPER IMPORTANT %%%
%%%%%%%%%%%%%%%%%%%%%%%
% TOFIX:
% - not entirely sure why all need to be multipled by -1. Was only
% expecting Vz...
% - something to do with x/y swap in MATLAB and that quiver often needs a
% dimension flip.

% Vx_masked = -1 .* Vx_masked;
% Vy_masked = -1 .* Vy_masked;
% Vz_masked = -1 .* Vz_masked;


%% Permute volumes to match MRtrix
% - I know vectors are accurate in MRtrix, so easiest to match with that

% % x = RL / y = SI --- NOT WORKING
% perm = [3,1,2,4];
% reori = @(vol) flip(flip(permute(vol,perm),1),2);
% Vx_masked = Vx_masked;
% Vy_masked = Vy_masked;
% Vz_masked = Vz_masked;

% x = RL / y = AP --- WORKING
perm = [2,1,3,4]; 
reori = @(vol) flip(flip(permute(vol,perm),1),2);
Vx_masked = -1 .* Vx_masked;
Vy_masked = -1 .* Vy_masked;
Vz_masked = -1 .* Vz_masked;

% % x = AP / y = SI --- NOT WORKING
% perm = [2,3,1,4]; 
% reori = @(vol) rot90(flip(permute(vol,perm),1));
% Vx_masked = Vx_masked;
% Vy_masked = Vy_masked;
% Vz_masked = Vz_masked;

CV = reori(cine_nii.img);
V.x = reori(Vx_masked);
V.y = reori(Vy_masked);
V.z = reori(Vz_masked);
V.m = reori(Vmag_masked);

% V.x = reori(-1.*Vz_masked);
% V.y = reori(-1.*Vx_masked);
% V.z = reori(-1.*Vy_masked);


%% View slice

sl = 35; %sl = sl + 1;   % +1 because MRtrix counts from 0
 t = 1;  %t = t + 1;    % +1 because MRtrix counts from 0

% meshgrid for quiver
[mx,my,mz] = meshgrid(1:size(V.x,2),1:size(V.x,1),1:size(V.x,3));

% quiver
% nb: imtar is 2D, quiver is 3D. Should probably use slice.m to match
% images and velocities on same grid.
imtar(CV(:,:,sl,t)); hold on;
q = quiver3( mx(:,:,1) , my(:,:,1), mz(:,:,1), ...
         V.x(:,:,sl,t), V.y(:,:,sl,t), V.z(:,:,sl,t) , ...
         'linewidth', 1, 'autoscalefactor', 2 );
quiverColour(q);
caxis([0 150]); colorbar('off'); axis('image'); % important that this comes last
title(['In MRtrix: Slice number = ' num2str(sl-1) ', Cardiac phase = ' num2str(t-1)]);
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);


%% Play slice as movie

% %TOFIX: colormap becomes gray - something to do with quiverColour I think
% %and also quiver3 changes the color each time you execute that line
% 
% nFrame = size(CV,4);
% 
% sl = 34; sl = sl + 1;   % +1 because MRtrix counts from 0
% 
% % meshgrid for quiver
% [mx,my,mz] = meshgrid(1:size(V.x,2),1:size(V.x,1),1:size(V.x,3));
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% for tt = 1:nFrame
% 
%     clf;
%     clear q
%     
%     t = tt;
% 
%     imagesc(CV(:,:,sl,t)); colormap(gray); axis image; hold on;
%     q = quiver3( mx(:,:,1) , my(:,:,1), mz(:,:,1), ...
%                  V.x(:,:,sl,t), V.y(:,:,sl,t), V.z(:,:,sl,t), ...
%                  'linewidth', 1, 'autoscalefactor', 2 );
%     quiverColour(q);
%     currentColormap = colormap(gca);
%     
%     title(['Frame number = ' num2str(tt)]);
% 
%     drawnow;
%     pause(0.5);
%     
% end


%% View single slice - slice.m version
% - true 3D, ie: vectors go through plane, but not great for visualisation

sl=35;
xslice = [];   
yslice = [];
zslice = [];

figure;
% q = quiver3( mx(:,:,sl)+0.5 , my(:,:,sl)+0.5, mz(:,:,sl)-1, ... % +0.5 to center on voxels. -3 so that vectors not hidden behind image.
%          V.x(:,:,sl,t), V.y(:,:,sl,t), V.z(:,:,sl,t) , ...
%          'linewidth', 1, 'autoscalefactor', 2 );
q = quiver( mx(:,:,sl)+0.5 , my(:,:,sl)+0.5, ... % +0.5 to center on voxels. -3 so that vectors not hidden behind image.
         V.x(:,:,sl,t), V.y(:,:,sl,t), ...
         'linewidth', 1, 'autoscalefactor', 2 );
quiverColour(q);

hold on;
% imslice = slice(CV(:,:,:,t), mx(:,:,sl) , my(:,:,sl), mz(:,:,sl), 'nearest'); % nearest interpolation matches imtar
imslice = slice(mx, my, mz, CV(:,:,:,t), xslice, yslice, zslice, 'nearest'); % nearest interpolation matches imtar

title(['In MRtrix: Slice number = ' num2str(sl-1) ', Cardiac phase = ' num2str(t-1)]);
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
% set(gca,'xticklabel',[]); set(gca,'xtick',[]); set(gca,'yticklabel',[]); set(gca,'ytick',[]);
colormap('gray'); grid('off'); set(imslice,'edgecolor','none'); axis('image');

view(2); camorbit(0,180); % same view as MRtrix



%% Generalised single slice - slice.m version
% - true 3D, ie: vectors go through plane, but not great for visualisation

% numVox.x = size(CV,1);
% numVox.y = size(CV,2);
% numVox.z = size(CV,3);
% 
% xslice = 35;   
% yslice = [];
% zslice = [];
% t = 1;
% 
% 
% % z-slice:
% if ~isempty(zslice)
%     xcoords = 1:numVox.x;   
%     ycoords = 1:numVox.y;
%     zcoords = zslice;
% 
%     figure;
% %     q = quiver3( mx(xcoords,ycoords,zcoords)+0.5 , my(xcoords,ycoords,zcoords)+0.5, mz(xcoords,ycoords,zcoords)-1, ... % +0.5 to center on voxels. -1 so that vectors not hidden behind image.
% %              V.x(xcoords,ycoords,zcoords,t), V.y(xcoords,ycoords,zcoords,t), V.z(xcoords,ycoords,zcoords,t) , ...
% %              'linewidth', 1, 'autoscalefactor', 2 );
%     q = quiver( mx(xcoords,ycoords,zcoords)+0.5 , my(xcoords,ycoords,zcoords)+0.5, ... % +0.5 to center on voxels. -3 so that vectors not hidden behind image.
%                 V.x(xcoords,ycoords,zcoords,t), V.y(xcoords,ycoords,zcoords,t), ...
%                 'linewidth', 1, 'autoscalefactor', 2 );
%     quiverColour(q);
% 
%     hold on;
%     % imslice = slice(CV(:,:,:,t), mx(:,:,sl) , my(:,:,sl), mz(:,:,sl), 'nearest'); % nearest interpolation matches imtar
%     imslice = slice(mx, my, mz, CV(:,:,:,t), xslice, yslice, zslice, 'nearest'); % nearest interpolation matches imtar
% 
%     title(['In MRtrix: Slice number = ' num2str(sl-1) ', Cardiac phase = ' num2str(t-1)]);
%     set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
%     set(gca,'xticklabel',[]); set(gca,'xtick',[]); set(gca,'yticklabel',[]); set(gca,'ytick',[]);
%     colormap('gray'); grid('off'); set(imslice,'edgecolor','none'); axis('image');
% 
%     view(2); camorbit(0,180); % same view as MRtrix
% 
% % x-slice --- NOT WORKING
% elseif ~isempty(xslice)
%     xcoords = xslice;   
%     ycoords = 1:numVox.y;
%     zcoords = 1:numVox.z;
% 
%     figure;
%     q = quiver3( mx(xcoords,ycoords,zcoords)+0.5 , my(xcoords,ycoords,zcoords)+0.5, mz(xcoords,ycoords,zcoords)-1, ... % +0.5 to center on voxels. -1 so that vectors not hidden behind image.
%              V.x(xcoords,ycoords,zcoords,t), V.y(xcoords,ycoords,zcoords,t), V.z(xcoords,ycoords,zcoords,t) , ...
%              'linewidth', 1, 'autoscalefactor', 2 );
% %     q = quiver( mx(xcoords,ycoords,zcoords)+0.5 , my(xcoords,ycoords,zcoords)+0.5, ... % +0.5 to center on voxels. -3 so that vectors not hidden behind image.
% %                 V.x(xcoords,ycoords,zcoords,t), V.y(xcoords,ycoords,zcoords,t), ...
% %                 'linewidth', 1, 'autoscalefactor', 2 );
%     quiverColour(q);
% 
%     hold on;
%     % imslice = slice(CV(:,:,:,t), mx(:,:,sl) , my(:,:,sl), mz(:,:,sl), 'nearest'); % nearest interpolation matches imtar
%     imslice = slice(mx, my, mz, CV(:,:,:,t), xslice, yslice, zslice, 'nearest'); % nearest interpolation matches imtar
% 
%     title(['In MRtrix: Slice number = ' num2str(sl-1) ', Cardiac phase = ' num2str(t-1)]);
%     set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
%     set(gca,'xticklabel',[]); set(gca,'xtick',[]); set(gca,'yticklabel',[]); set(gca,'ytick',[]);
%     colormap('gray'); grid('off'); set(imslice,'edgecolor','none'); axis('image');
% 
%     view(2); camorbit(0,180); % same view as MRtrix
% end
