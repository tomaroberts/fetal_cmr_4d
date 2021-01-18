%% Calculate Coil Sensitivity Maps from Dynamic Undersampled kt data

%% Admin
reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr254_tom_recon';
ktDir = '\ktrecon';

stackNum = 20;
ktFactor = 8;
sliceNum = 5;
shotNum  = 1; % 1:nDyn/ktFactor


%% Anon Fns
dimX = 1;
dimY = 2;

kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, dimX ), dimY ) );

phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [+1;-1], size(k,1)/2, 1 ) ) ) );
 
inv_phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [-1;+1], size(k,1)/2, 1 ) ) ) );

fn_rmv_os = @( x ) x((round(size(x,dimX)/4)+1):(size(x,dimX)-round(size(x,dimX)/4)),:,:,:);


%% Load data
cd([reconDir ktDir]);
load(['s' num2str(stackNum) '_csm.mat']);
load(['s' num2str(stackNum) '_kspace.mat'],'ktAcq');

% Load unprocessed csm data
cd([reconDir, '/csm']); %nb: only for 254 atm
C0 = load('csmFile_SENS_0_except_CalcSensitivity.mat');
csm0 = C0.csm; imBody0 = C0.imBody; imCoil0 = C0.imCoil; 
clear C0


%% View mRecon-derived Coil Sensitivity Maps
montage_RR(squeeze(csm(:,:,1,:,sliceNum)),['Coil Sens. Maps, slice' num2str(sliceNum)],[0,1]); colormap('jet');
montage_RR(abs(imCoil(:,:,:,:,sliceNum)),['Array Coils, slice' num2str(sliceNum)],[0 5e4]);
montage_RR(abs(imBody(:,:,:,:,1)),['Body Coil, slice' num2str(sliceNum)],[0 5e4]);
 
% %% Sanity check: 
% % csm = imCoil ./ imBody
% for tt = 1:size(imCoil,4)
%     csmTom(:,:,1,tt,:) = imCoil(:,:,1,tt,:) ./ imBody(:,:,1,:,:);
% end
% montage_RR(squeeze(csmTom(:,:,1,:,sliceNum)),['Coil Sens. Maps, slice' num2str(sliceNum)],[0,1]); colormap('jet');
% csmTom is basically = csm, barring some masking and tidying in MRecon


%% Create Single Shot k-spaces
% Use these in place of imCoil?

dim.x    = size(ktAcq,1);
dim.y    = size(ktAcq,2);
dim.dyn  = size(ktAcq,3);
dim.chan = size(ktAcq,4);
dim.loca = size(ktAcq,5);

ktAcqSS = zeros(dim.x, dim.y, ktFactor, dim.dyn/ktFactor, dim.chan, dim.loca);
ktAcqSS = reshape(ktAcq , dim.x, dim.y, ktFactor, dim.dyn/ktFactor, dim.chan, dim.loca );

ktAcqSS = phase_correct( ktAcqSS );

xtAcqSS = kt2xt( ktAcqSS );

% memory clearance
clear ktAcq

% View Single Shot Data from All Coils
montage_RR(abs(squeeze(sum(xtAcqSS(:,:,:,shotNum,:,sliceNum),3))),'',[0,1]); title(['Superframe number ' num2str(shotNum)]);
montage_RR(abs(imCoil0(:,:,:,:,sliceNum)),'',[0 5e4]); title('No MRecon Processing - imCoil');


%% Truncated k-space
% i.e.: synthetic training data
trcWidth = 11;
midK     = dim.y/2;
trcRange = midK - floor(trcWidth/2) : midK + floor(trcWidth/2);

ktAcqTrc = zeros(dim.x, dim.y, dim.dyn/ktFactor, dim.chan, dim.loca);
ktAcqTrc(:,trcRange,:,:,:) = squeeze(sum(ktAcqSS(:,trcRange,:,:,:,:),3));

xtAcqTrc = kt2xt( ktAcqTrc );

% View
imtar(abs(squeeze(ktAcqTrc(:,:,shotNum,1,sliceNum))),0,100); % just one coil
montage_RR(abs(squeeze(xtAcqTrc(:,:,shotNum,:,sliceNum))),'',[0,1]); title(['Central ' num2str(trcWidth) ' lines of k-space']);

% Average across multiple Superframes
sfRange = 1:12;
montage_RR(abs(squeeze(mean(xtAcqTrc(:,:,sfRange,:,sliceNum),3))),'',[0,1]); title(['Central ' num2str(trcWidth) ' lines of k-space // Averaged over Superframes: ' num2str(sfRange)]);



%% SoS on images
% Use in place of imBody?
xtAcqShot = zeros(dim.x, dim.y, dim.dyn/ktFactor, 1, dim.loca);
for iShot = 1:(dim.dyn/ktFactor)
    for iSlice = 1:dim.loca
        xtAcqSoS(:,:,iShot,:,iSlice) = sqrt( sum( abs( sum(xtAcqSS(:,:,:,iShot,:,iSlice),3) ).^2 ,5) );
        xtAcqTrcSoS(:,:,iShot,:,iSlice) = sqrt( sum( abs( xtAcqTrc(:,:,iShot,:,iSlice) ).^2 ,4) );
    end
end
montage_RR(squeeze(xtAcqSoS(:,:,shotNum,:,sliceNum)),'',[0,5]);
montage_RR(squeeze(xtAcqTrcSoS(:,:,shotNum,:,sliceNum)),'',[0,5]);

% DC SoS Image
for iSlice = 1:dim.loca
    xtDcSoS(:,:,iSlice) = squeeze( sum( xtAcqSoS(:,:,:,:,iSlice), 3));
end
montage_RR(xtDcSoS(:,:,sliceNum),'',[0,50] );

% Body mask
BM = ~( xtDcSoS < prctile( xtDcSoS(:), 75 ) );
for iSlice = 1:dim.loca
    BMfilt(:,:,iSlice) = logical( medfilt2( double( BM(:,:,iSlice) ))); % remove speckle
end
imtar([BM(:,:,sliceNum) BMfilt(:,:,sliceNum)]); 
BM = BMfilt; clear BMfilt



%% Jo threshold
idx = xtAcqSoS > 0.05;
xtAcqSoSThr = xtAcqSoS .* idx;



%% Calculate Shot-by-shot Coil Sensitivity Maps
csmSS = zeros(size(xtAcqSS,1),size(xtAcqSS,2),size(xtAcqSS,4),1,size(xtAcqSS,6),size(xtAcqSS,5));
csmSSsmooth = zeros(size(csmSS));

csmTrc = zeros(dim.x, dim.y, dim.dyn/ktFactor, 1, dim.loca);

for iShot = 1:(dim.dyn/ktFactor)
    for iSlice = 1:dim.loca
        for iCoil = 1:dim.chan
            
            % csmSS = xtAcqSS ./ xtAcqSoS
            % i.e.: csmSS = virtual imCoil ./ virtual imBody
%             csmSS(:,:,iShot,:,iSlice,iCoil) = BM(:,:,iSlice) .* ... % Apply mask            
            csmSS(:,:,iShot,:,iSlice,iCoil) = ... % No mask
                squeeze( abs( sum( xtAcqSS( :,:,:,iShot,iCoil,iSlice ), 3) ) ) ./ xtAcqSoS(:,:,iShot,:,iSlice) ;
            
            % smooth with median filter
            csmSSsmooth(:,:,iShot,:,iSlice,iCoil) = medfilt2(csmSS(:,:,iShot,:,iSlice,iCoil) , [5,5]);
        
            % Truncated version
            csmTrc(:,:,iShot,:,iSlice,iCoil) = BM(:,:,iSlice) .* ...
                squeeze( abs( sum( xtAcqTrc( :,:,iShot,iCoil,iSlice ), 3) ) ) ./ xtAcqTrcSoS(:,:,iShot,:,iSlice);
 
            
            
            
        end
    end
end

cMax = 0.25; % nb: need to figure out why doesn't look good on 0-1 scale, as for original csm
montage_RR(squeeze(csmSS(:,:,shotNum,:,sliceNum,:)),['Single-shot Coil Sens. Maps, shot' num2str(shotNum)],[0,cMax]); colormap('jet');
montage_RR(squeeze(csmSSsmooth(:,:,shotNum,:,sliceNum,:)),['Single-shot Coil Sens. Maps, shot' num2str(shotNum)],[0,cMax]); colormap('jet');

montage_RR(squeeze( sum( csmSS(:,:,:,:,sliceNum,:), 3) ),['Sum of Single-shot Coil Sens. Maps'],[0,cMax*1e1]); colormap('jet');

% Unprocessed imCoil vs. truncated undersampled k-space
montage_RR(squeeze( sum( csm0(:,:,:,:,sliceNum,:), 3) ),'',[0,1]); colormap('jet'); title('Unprocessed scanner csm');
montage_RR(squeeze( sum( csmTrc(:,:,:,:,sliceNum,:), 3) ),'',[0,3]); colormap('jet'); title(['Truncated synthetic csm // ' num2str(trcWidth) ' lines of k-space']);


%% Gauss filtered csm
% - might be unnecessary as added resolution beneficial for csm?
% - mainly testing to see if can get looking like original csm
filtSz = 8;
csmSSfilt = zeros(size(csmSS));

for iShot = 1:(dim.dyn/ktFactor)
    for iSlice = 1:dim.loca
        for iCoil = 1:dim.chan
            csmSSfilt(:,:,iShot,:,iSlice,iCoil) = imgaussfilt( csmSS(:,:,iShot,:,iSlice,iCoil), filtSz );
        end
    end
end

montage_RR(squeeze(csmSSfilt(:,:,shotNum,:,sliceNum,:)),['Filtered single-shot Coil Sens. Maps, shot' num2str(shotNum)],[0,0.25]); colormap('jet');


