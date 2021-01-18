% fCohort = [189 191 194 197 201 202 214];

fCohort = [194 197 201 202 214];

for ff = 1:5
    
    fNum = fCohort(ff);
    
    disp(['Running fcmr' num2str(fNum) ' ... ' ]);
    
    fcmrDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data';
    
    cd([fcmrDir '\fcmr' num2str(fNum)]);
    
    velDir = 'vel_vol';
    cd(velDir);
    
    
    
    %% load vel
    if ~isfile('velocity-final-RESLICE-0.nii.gz') || ~isfile('velocity-final-RESLICE-1.nii.gz') || ~isfile('velocity-final-RESLICE-2.nii.gz')
        reslice_nii('velocity-final-0.nii.gz', 'velocity-final-RESLICE-0.nii.gz');
        reslice_nii('velocity-final-1.nii.gz', 'velocity-final-RESLICE-1.nii.gz');
        reslice_nii('velocity-final-2.nii.gz', 'velocity-final-RESLICE-2.nii.gz');
    end
    
    velx_nii = load_nii('velocity-final-RESLICE-0.nii.gz');
    vely_nii = load_nii('velocity-final-RESLICE-1.nii.gz');
    velz_nii = load_nii('velocity-final-RESLICE-2.nii.gz');
    
    velx_nii.img = velx_nii.img .* 1e2;
    vely_nii.img = vely_nii.img .* 1e2;
    velz_nii.img = velz_nii.img .* 1e2;
    
    cd([fcmrDir '\fcmr' num2str(fNum)]);
    cd(velDir);
    mkdir vel_drift_comparison
    cd vel_drift_comparison
    
    % save 4D
    Vx3D_nii = load_untouch_nii([fcmrDir '\fcmr' num2str(fNum) '\' velDir '\velocity-final-RESLICE-0.nii.gz']);
    
    for tt = 1:25
        v3D = Vx3D_nii;
        v3D.img = cat(4,velx_nii.img(:,:,:,tt), vely_nii.img(:,:,:,tt), -1 .* velz_nii.img(:,:,:,tt));
        
        v3D.hdr.dime.dim(1) = 4;
        v3D.hdr.dime.dim(2:5) = [size(velx_nii.img,1) size(velx_nii.img,2) size(velx_nii.img,3) 3];
        v3D.hdr.dime.pixdim(5) = 0;
        
        save_untouch_nii(v3D,['vel_vol_noCorr' num2str(tt-1) '.nii.gz']);
    end
    
    disp('Finished making noCorr .nii files ...');
    
    
    
    %% load polyCorr_vell
    cd([fcmrDir '\fcmr' num2str(fNum)]);
    cd(velDir);
    
    velx_nii = load_nii('velocity-final-polyCorr-RESLICE-0.nii.gz');
    vely_nii = load_nii('velocity-final-polyCorr-RESLICE-1.nii.gz');
    velz_nii = load_nii('velocity-final-polyCorr-RESLICE-2.nii.gz');
    
    velx_nii.img = velx_nii.img .* 1e2;
    vely_nii.img = vely_nii.img .* 1e2;
    velz_nii.img = velz_nii.img .* 1e2;
    
    cd([fcmrDir '\fcmr' num2str(fNum)]);
    cd(velDir);
    mkdir vel_drift_comparison
    cd vel_drift_comparison
    
    % save 4D
    Vx3D_nii = load_untouch_nii([fcmrDir '\fcmr' num2str(fNum) '\' velDir '\velocity-final-RESLICE-0.nii.gz']);
    
    for tt = 1:25
        v3D = Vx3D_nii;
        v3D.img = cat(4,velx_nii.img(:,:,:,tt), vely_nii.img(:,:,:,tt), -1 .* velz_nii.img(:,:,:,tt));
        
        v3D.hdr.dime.dim(1) = 4;
        v3D.hdr.dime.dim(2:5) = [size(velx_nii.img,1) size(velx_nii.img,2) size(velx_nii.img,3) 3];
        v3D.hdr.dime.pixdim(5) = 0;
        
        save_untouch_nii(v3D,['vel_vol_polyCorr' num2str(tt-1) '.nii.gz']);
    end
    
    disp('Finished making polyCorr .nii files ...');
    
end