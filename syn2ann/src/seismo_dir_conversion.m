function [varargout] = seismo_dir_conversion(varargin)
    cpp = varargin{1};
    
    if any(strcmpi(cpp,{'e';'w';'ew';'we';'wec'}))
        cpn = {'ew'};
    end
    
    if any(strcmpi(cpp,{'n';'s';'ns';'sn';'nsc'}))
        cpn = {'ns'};
    end
    
    if any(strcmpi(cpp,'gh'))
        cpn = {'ew';'ns'};
    end
    
    if any(strcmpi(cpp,{'z';'ud';'udc'}))
        cpn = {'ud'};
    end
    
    varargout{1} = cpn;
    return
end
