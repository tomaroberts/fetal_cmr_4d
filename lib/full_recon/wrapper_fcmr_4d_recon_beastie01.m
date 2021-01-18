
scanDate = '2017_09_07';
patID    = 'WI_374747';
reconDir = '/scratch/cmo19/Data/fcmr202';

cd(reconDir);

fcmr_4d_recon_beastie01('scanDate',scanDate,'patID',patID,'reconDir',reconDir);


%% c_fcmr375

scanDate   = '2020_12_11';
patID      = 'EL_421313';
reconDir   = '/scratch/cmo19/Data/c_fcmr375';
dlgChoice  = 4;
rawDataDir = '/pnraw01/raw-ingenia'; % or archive
seriesNos  = [20 21 22 23 24];

fcmr_4d_recon_beastie01('scanDate',scanDate, 'patID',patID, 'reconDir',reconDir, 'dlgChoice', dlgChoice);





%% Without specifying

fcmr_4d_recon_beastie01;

