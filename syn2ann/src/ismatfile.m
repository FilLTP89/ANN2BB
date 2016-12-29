function ismat=ismatfile(str)
    
    try
        ismat=~isempty(whos('-file',str));
    catch
        try 
            ismat=~isempty(whos('-file',str(1:end-4)));
        catch
            ismat=false;
        end
    end
    
    return
end