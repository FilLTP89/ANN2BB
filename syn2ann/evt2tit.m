function [varargout] = evt2tit(varargin)
    evt = varargin{1};
    typ = varargin{2};
    
    switch lower(typ)
        case 'itaca'
            evt(strfind(evt,'.')) = '';
            evt = evt(1:end-2);
            fgn = evt(3:end);
        case {'kiknet','knet'}
            fgn = evt;
            evt = strcat('20',evt);
    end
    evt = strcat(evt(1:4),'-',evt(5:6),'-',evt(7:8),{' '},evt(9:10),':',evt(11:12));
    
    varargout{1} = evt;
    varargout{2} = fgn;
    return
end