%% Cribbed from: fcmr_4dflow_vessel_analysis

function ROI = fcmr_4dflow_roi2vel(ROI, Vsum, Vmag, nFrame)

for tt = 1:nFrame
    ROI.vsum(tt) = mean(nonzeros(Vsum(:,:,:,tt) .* ROI.img(:,:,:,tt)));
    ROI.vmag(tt) = mean(nonzeros(Vmag(:,:,:,tt) .* ROI.img(:,:,:,tt)));
    ROI.vstd(tt) = std(nonzeros(Vmag(:,:,:,tt) .* ROI.img(:,:,:,tt)));
    ROI.vmax(tt) = max(nonzeros(Vmag(:,:,:,tt) .* ROI.img(:,:,:,tt)));
    ROI.vmin(tt) = min(nonzeros(Vmag(:,:,:,tt) .* ROI.img(:,:,:,tt)));
end

% rectify vsum
ROI.vsum = sign( sum( ROI.vsum(:) ) ) * ROI.vsum;

end