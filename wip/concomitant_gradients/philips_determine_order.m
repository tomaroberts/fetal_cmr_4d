function order = philips_determine_order(str_XYZ)

if strcmp(str_XYZ,'+x')
    order = +1;
elseif strcmp(str_XYZ,'-x')
    order = -1;
elseif strcmp(str_XYZ,'+y')
    order = +2;
elseif strcmp(str_XYZ,'-y')
    order = -2;
elseif strcmp(str_XYZ,'+z')
    order = +3;
elseif strcmp(str_XYZ,'-z')
    order = -3;
end



end