function ismat=ismatfile(str)
    
    try
        ismat=~isempty(whos('-file',str));
    catch
        ismat=false;
    end
    
    return
end