% qflow
qDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow\data_ffd';
cd(qDir);

velFiles = dir('*_vel_*scaled.nii.gz');
velFiles = {velFiles([1,2,3]).name}; %re-order: ap/fh/rl

for ii = 1:numel(velFiles)    
    nii = load_nii(velFiles{ii});   
    VEL(:,:,:,ii) = nii.img;
end

vq = double(VEL);

% mask
M = load_nii('../mask/inside_pipes_mask.nii.gz');
M = single(M.img);

implay_RR([M.*vq(:,:,:,1), ...
           M.*vq(:,:,:,2), ...
           M.*vq(:,:,:,3)] ...
           ,'jet',[-5,5]);


% bssfp
bDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\data_ffd';
cd(bDir);

nii = load_nii('ap_vel_ffd.nii.gz');
vb(:,:,:,1) = double(nii.img);
nii = load_nii('fh_vel_ffd.nii.gz');
vb(:,:,:,2) = -1.*double(nii.img);
nii = load_nii('rl_vel_ffd.nii.gz');
vb(:,:,:,3) = double(nii.img);


% masks
cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\mask');
P = load_nii('inside_pipes_mask.nii.gz');
S = load_nii('spherical_part_mask.nii.gz');

% W = load_nii('static_water_mask_128px_ffd.nii.gz');
% B = load_nii('background_mask.nii_128px.gz');
% B = load_nii('background_and_artefact_mask_128px.nii.gz');
% M = double(P.img & ~W.img & ~B.img);

M = double(P.img & S.img);

vq3d = sqrt(vq(:,:,:,1).^2 + vq(:,:,:,2).^2 + vq(:,:,:,3).^2);
vb3d = sqrt(vb(:,:,:,1).^2 + vb(:,:,:,2).^2 + vb(:,:,:,3).^2);

implay_RR([M.*vq(:,:,:,1), M.*vq(:,:,:,2), M.*vq(:,:,:,3), M.*vq3d;...
           M.*vb(:,:,:,1), M.*vb(:,:,:,2), M.*vb(:,:,:,3), M.*vb3d]...
           ,'jet',[-10,10]);

       
       

%% Plot histograms

lowerThresh = -20;
upperThresh = 20;
binSize = 0.2;

bkrdClip = 0.01;

for ii = 1:3
    M_vq(:,:,:,ii) = M.*vq(:,:,:,ii);
    M_vb(:,:,:,ii) = M.*vb(:,:,:,ii);
end

M_vq3d = M.*vq3d;
M_vb3d = M.*vb3d;

figure;

% set background/really low values to NaN
idx = find(M_vq>-bkrdClip & M_vq<bkrdClip); M_vq(idx) = NaN;
idx = find(M_vb>-bkrdClip & M_vb<bkrdClip); M_vb(idx) = NaN;
idx = find(M_vq3d>-bkrdClip & M_vq3d<bkrdClip); M_vq3d(idx) = NaN;
idx = find(M_vb3d>-bkrdClip & M_vb3d<bkrdClip); M_vb3d(idx) = NaN;

% % display with thresholding
% implay_RR([M_vq(:,:,:,1), M_vq(:,:,:,2), M_vq(:,:,:,3), M_vq3d;...
%            M_vb(:,:,:,1), M_vb(:,:,:,2), M_vb(:,:,:,3), M_vb3d]...
%            ,'jet',[-10,10]);


for ii = 1:3
    subplot(2,2,ii);
    temp = M_vq(:,:,:,ii); disp(num2str(numel(temp)));
    QVhist = temp(:); %QVhist = QVhist(QVhist>lowerThresh); QVhist = QVhist(QVhist<upperThresh);
    histogram(QVhist,'binEdges',lowerThresh:binSize:upperThresh,'Normalization','Probability');

    hold on;

    temp = M_vb(:,:,:,ii); disp(num2str(numel(temp)));
    BVhist = temp(:); %QVhist = QVhist(QVhist>lowerThresh); QVhist = QVhist(QVhist<upperThresh);
    histogram(BVhist,'binEdges',lowerThresh:binSize:upperThresh,'Normalization','Probability');

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
end

lowerThresh = 0.5;
upperThresh = 20;
binSize = 0.1;

subplot(2,2,4);
QVhist = M_vq3d(:); %QVhist = QVhist(QVhist>lowerThresh); QVhist = QVhist(QVhist<upperThresh);
histogram(QVhist,'binEdges',0:binSize:upperThresh,'Normalization','Probability');

hold on;

BVhist = M_vb3d(:); %QBhist = QBhist(QBhist>lowerThresh); QBhist = QBhist(QBhist<upperThresh);
histogram(BVhist,'binEdges',0:binSize:upperThresh,'Normalization','Probability');

legend('PC-SPGR','PC-bSSFP');
xlabel('Velocity [cm/s]');
ylabel('Counts');
title('Magnitude');

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


