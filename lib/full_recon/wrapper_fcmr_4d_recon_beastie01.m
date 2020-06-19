
scanDate = '2017_09_07';
patID    = 'WI_374747';
reconDir = '/scratch/cmo19/Data/fcmr202';

cd(reconDir);

fcmr_4d_recon_beastie01('scanDate',scanDate,'patID',patID,'reconDir',reconDir);


%% c_fcmr273

scanDate = '2019_08_23';
patID    = 'RI_128958';
reconDir = '/scratch/cmo19/Data/c_fcmr273';

fcmr_4d_recon_beastie01('scanDate',scanDate,'patID',patID,'reconDir',reconDir);





%% Without specifying

fcmr_4d_recon_beastie01;

