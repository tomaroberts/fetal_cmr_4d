function str_XYZ = philips_lpf2xyz(str_LPF)

% converts from Patient System (LPF = LR-PA-FH) 
% to XYZ system on Philips scanner

if strcmp(str_LPF,'PA')
    str_XYZ = '+x';
elseif strcmp(str_LPF,'AP')
    str_XYZ = '-x';
elseif strcmp(str_LPF,'RL')
    str_XYZ = '+y';
elseif strcmp(str_LPF,'LR')
    str_XYZ = '-y';
elseif strcmp(str_LPF,'FH')
    str_XYZ = '+z';
elseif strcmp(str_LPF,'HF')
    str_XYZ = '-z';
end

end