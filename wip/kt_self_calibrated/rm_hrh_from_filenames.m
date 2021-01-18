
% Remove _hrh from filenames

seriesNos = [20 21 22 24];

cd(reconDir);
for seriesNo = seriesNos

    x = dir(['ktrecon_hrh/s' num2str(seriesNo) '*hrh*']);
    x(5)=[]; x(5)=[]; %ignores _hrh_recon.mat and _kspace_hrh.mat (Lucilio outputs)
    
    cd ktrecon_hrh
    for ii = 1:numel(x)
        movefile(x(ii).name,erase(x(ii).name,'_hrh'))
    end
    cd ..
    
end




% seriesNos = [14 15 16 17 18];
% 
% cd(reconDir);
% cd ktrecon
% for seriesNo = seriesNos
% 
% %     load(['s' num2str(seriesNo) '_dc_recon.mat'],'xtHrhDc');
% %     xtDc = xtHrhDc;
% %     save(['s' num2str(seriesNo) '_dc_recon.mat'],'xtDc');
%     
%     load(['s' num2str(seriesNo) '_rlt_recon.mat'],'xtHrh');
%     xtRcn = xtHrh;
%     save(['s' num2str(seriesNo) '_rlt_recon.mat'],'xtRcn');
%     
% end