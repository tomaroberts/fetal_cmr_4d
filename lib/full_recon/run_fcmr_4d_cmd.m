
%% Run 4D FCMR pipeline

scanDate   = '2020_12_04';
patID      = 'KI_415285';
reconDir   = '/scratch/cmo19/Data/c_fcmrfimox028air';
dlgChoice  = 4;
rawDataDir = '/pnraw01/raw-ingenia'; % or archive
seriesNos  = [13 14 15 16 17];

fcmr_4d_recon_beastie01('scanDate',  scanDate,   ...
                        'patID',     patID,      ...
                        'reconDir',   reconDir,   ...
                        'dlgChoice',  dlgChoice,  ...
                        'rawDataDir', rawDataDir, ...
                        'seriesNos',  seriesNos   ...
                       );
