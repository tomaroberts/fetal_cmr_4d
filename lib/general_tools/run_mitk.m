function run_mitk( )

%% RUN_MITK
%
% Quick function to call MITK Workbench using system command
%
% For use on beastie01
%
% Note: MITK can't be called from within Matlab 2015a on beastie01
% Hence, this _horrible_ work around calls Matlab 2016b within the
% run_mitk.bash script.
%
%

% Path to MITK Workbench bash script
mitkWkbhScriptPath = '/home/cmo19/MATLAB/run_mitk.bash';

% Load MITK Workbench
fprintf('\nLoading MITK Workbench ...\n');
cmdStr = sprintf('''%s''', mitkWkbhScriptPath );
[~,~] = system(cmdStr);
fprintf('Exiting MITK Workbench.\n');

% end run_mitk(...)
end