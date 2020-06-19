%% 2019_07_04 --- PC-bSSFP --- FlowPhantom_9 --- spherical flask
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1)
%% Files reconstructed on beastie02 using FP9_pc_mrecon_from_raw.m
%- outputs: _RE / _IM / _CPLX / _MAG / _PH.nii.gz files + others
reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP';
cd(reconDir);

%% Step 2)
%% Make directory stucture and copy magnitude/complex data
% In bash:
% mkdir data
% mkdir mask
%
% cp *_MAG* ../data/
% cp *_CPLX* ../data/
% cp *_RE* ../data/
% cp *_IM* ../data/
%
% Rename accordingly...
%
% Make ../mask/static_water_mask.nii.gz in MITK
%
%
% AT THIS POINT, I resampled the bssfp data to 128x128x120 px to match the
% QFlow data. Script for this:
% ../data_128px/bssfp_downsample.bash

%% Step 3)
%% Make phase data
% cd(reconDir);
% cd data
% fileNames = dir('*_re.nii.gz');
% 
% for ii = 1:numel(fileNames)
%     currStr = erase(fileNames(ii).name, '_re.nii.gz');
%     
%     re = load_untouch_nii( [currStr '_re.nii.gz'] );
%     im = load_untouch_nii( [currStr '_im.nii.gz'] );
%     
%     % initialise ph/ph_corr nifti structures
%     ph = re; ph_corr = re;
%         
%     % Create uncorrected phase images
%     cx.img = re.img + 1i.*im.img;
%     ph.img = angle(abs(cx.img).*exp(sqrt(-1)*(angle(cx.img)+pi))); %DO I NEED +pi in exp? Not sure
% 
%     % save ph (without polynomial correction)
%     ph.img = single(ph.img); % convert to single to match original .nii
%     save_untouch_nii(ph,[currStr '_ph.nii.gz']);
% 
%     % Create polynomial corrected phase images
%     cd ../mask
% %     static_water_mask = load_untouch_nii('static_water_mask.nii.gz');
%     static_water_mask = load_untouch_nii('qflow_static_water_mask.nii.gz');
%     static_water_mask.img = imgaussfilt3(static_water_mask.img,2.5);
%     
%     % run polynomial correction
%     % - Y is corrected complex data, ie: cx_corr.img
%     [Y, P0, P1] = phase_correction_poly_phantom( cx.img, ...
%                                'staticwatermask', logical(static_water_mask.img) );
%     ph_corr.img = angle(abs(Y).*exp(sqrt(-1)*(angle(Y))));
% 
%     % save ph_corr.nii.gz
%     ph_corr.img = single(ph_corr.img); % convert to single to match original .nii
%     
%     cd ../data
%     save_untouch_nii(ph_corr,[currStr '_ph_corr.nii.gz']);
% end

%% Step 4)
%% Register stacks into same space
% FFD: ../reconDir/data_ffd/register_bssfp.bash

%% Step 5)
%% Load stacks
reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP';
cd(reconDir);
cd data_ffd

magFiles = dir('*_mag_*.nii.gz');
magFiles = {magFiles([1,3,5,2,4,6]).name}; %re-order: ax/cor/sag/axrot/corrot/sagrot

phFiles = dir('*_ph_corr_*.nii.gz');
phFiles = {phFiles([1,3,5,2,4,6]).name};

for ii = 1:numel(magFiles)
    nii = load_untouch_nii(magFiles{ii});
    MAG(:,:,:,ii) = nii.img;
    
    nii = load_untouch_nii(phFiles{ii});
    PH(:,:,:,ii) = nii.img;
end

clear nii

% implay_RR([MAG(:,:,:,1), MAG(:,:,:,2), MAG(:,:,:,3);...
%            MAG(:,:,:,4), MAG(:,:,:,5), MAG(:,:,:,6)]...
%            ,'gray',[0,10]);
% 
% implay_RR([PH(:,:,:,1), PH(:,:,:,2), PH(:,:,:,3);...
%            PH(:,:,:,4), PH(:,:,:,5), PH(:,:,:,6)]...
%            ,'gray',[-pi/6,pi/6]);

%% Step 6)
%% Get gradient first moments
cd(reconDir);
cd raw

% first moments from GVE
Vm = 16.43;
Vp = 0;
Vs = -36.63;
Vmps = [Vm; Vp; Vs];

goalcFiles = dir('*_goalc.txt'); %already correct order: ax/cor/sag/axrot/corrot/sagrot

disp('Getting gradint first moments ...');
for ii = 1:numel(goalcFiles)
    gc = get_pcmr_orientation_parameters( goalcFiles(ii).name );
    [Vworld(:,ii), Vxyz(:,ii)] = gradfirstmom_mps2world(Vmps,gc);
end
disp('Got moments.');

cd ../data_ffd

%% Step 7)
%% Solve for Velocity
clear P Pmat V Vxyz_vec Vworld_vec VEL

% Choose stacks to use in reconstruction
% 1 = ax      /  2 = cor      /  3 = sag
% 4 = ax_rot  /  5 = cor_rot  /  6 = sag_rot
stackChoices = [1:6];

for ii = 1:size(PH,4)
    temp = PH(:,:,:,ii);
    PH_vec(:,ii) = temp(:);
end; clear temp

gamma = 2 .* pi .* 42577; %Hz/mT
Mscaling = (1e-3).^2;  %ms^2.mT/m --- First Moment scaling into correct units

V_xyz   = gamma .* Mscaling .* Vxyz(:,stackChoices)';
V_world = gamma .* Mscaling .* Vworld(:,stackChoices)';
Pmat    = PH_vec(:,stackChoices)';

% Solve
disp('Solving ...');
for ii = 1:length(Pmat)
    Vxyz_vec(:,ii) = V_xyz\Pmat(:,ii); %X = A\B => Vxyz*Vxyz_vec = Pmat
    Vworld_vec(:,ii) = V_world\Pmat(:,ii);
end
disp('Solved.');


%% Step 8)
%% Make Velocity Volume

% Un-vectorise and convert to velocity stacks
subV=ind2subV([size(PH,1),size(PH,2),size(PH,3)],1:length(Pmat));

disp('Making velocity volumes ...');
for ii = 1:length(subV)
    i = subV(ii,1);
    j = subV(ii,2);
    k = subV(ii,3);
    
    %xyz coords
    VEL.vx(i,j,k) = Vxyz_vec(1,ii) .* 1e2; %cm/s
    VEL.vy(i,j,k) = Vxyz_vec(2,ii) .* 1e2;
    VEL.vz(i,j,k) = Vxyz_vec(3,ii) .* 1e2;
    
    %world coords
    VEL.rl(i,j,k) = Vworld_vec(1,ii) .* 1e2;
    VEL.ap(i,j,k) = Vworld_vec(2,ii) .* 1e2;
    VEL.fh(i,j,k) = Vworld_vec(3,ii) .* 1e2;
end
disp('Created velocity volumes.');


%% View

M = load_untouch_nii('../mask/inside_pipes_mask.nii.gz');
M = double(M.img);

permOrder = [2,1,3];

implay_RR(permute(M.*MAG(:,:,:,1),permOrder),'gray');

implay_RR([M.*VEL.ap, M.*VEL.fh, M.*VEL.rl ...
           ],'jet',[-10,10]);
       
save('../FP9_bSSFP.mat','VEL','M');
      

%% save velocity NII

cd(reconDir);
cd data_ffd

templateNii = 'ax_mag_ffd.nii.gz';
nii = load_untouch_nii(templateNii);
nii.hdr.dime.glmax = 30;
nii.hdr.dime.glmin = -30;

mkdir('6_phase_stacks_recon');
cd('6_phase_stacks_recon');

nii.img = VEL.rl;
save_untouch_nii(nii,'rl_vel_ffd.nii.gz');
nii.img = VEL.ap;
save_untouch_nii(nii,'ap_vel_ffd.nii.gz');
nii.img = -1 .* VEL.fh; %%% IMPORTANT FOR MRtrix - must flip z/fh direction
save_untouch_nii(nii,'fh_vel_ffd.nii.gz');


