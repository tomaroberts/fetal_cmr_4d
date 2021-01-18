
cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr189_hrh\ktrecon');

load('s22_hrh_recon.mat');
xtTrn = squeeze(xtRcn);

dimX = 1;
fn_rmv_os = @( x ) x((round(size(x,dimX)/4)+1):(size(x,dimX)-round(size(x,dimX)/4)),:,:,:);

dimT = 3;
xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );
xf2xt = @( xf ) ifft( ifftshift( xf, dimT ), [], dimT );

xtTrn = fn_rmv_os(xtTrn);

sl = 7;
xtTrn = xtTrn(:,:,:,sl);
xtTrn = xtTrn./max(xtTrn(:)); % dirty normalise
xfTrn = xt2xf(xtTrn);
xfTrn(:,:,49) = 0;

imtar(xtTrn(:,:,1) ,0,.25)

imtar(permute(xfTrn(:,53:2:61,:),[1,3,2]),0,1); colormap('parula');


%% View xf

xfLine = 57;

figure; plot(squeeze(abs(xfTrn(95,xfLine,:))),'-o'); axis([0 96 0 0.25])



%% Find DC aliases

% dcAliases = 1:12:96;
aPadding  = 3;
dcAliases = 1:12:96+aPadding;
dcMin = dcAliases-aPadding;
dcMax = dcAliases+aPadding;

dcMin(dcMin<1)=1;
dcMax(dcMax>96)=96;

dcMask = zeros(200,96);
for dd = 1:(8+1)
    dcMask(:,dcMin(dd):dcMax(dd)) = 1;
end
   

imtar(dcMask .* ...
    permute(xfTrn(:,xfLine,:),[1,3,2]) ...
    ,0,1); colormap('parula');

dcAliasSigs = mean( dcMask .* permute(xfTrn(:,xfLine,:),[1,3,2]) ,1 );

figure; plot(abs(dcAliasSigs));


