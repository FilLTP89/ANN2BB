function [varargout] = evt2tit(varargin)
    evt = varargin{1};
    typ = varargin{2};
    
    switch lower(typ)
        case 'itaca'
            evt(strfind(evt,'.')) = '';
            fgn = evt(1:end);
            evt = evt(1:end-2);
        case {'kiknet','knet'}
            fgn = evt;
            evt = strcat('20',evt);
    end
    evt = strcat(evt(1:4),'-',evt(5:6),'-',evt(7:8),{' '},evt(9:10),':',evt(11:12));
    fgn = strcat(fgn(1:8),'_',fgn(9:14));
    varargout{1} = evt;
    varargout{2} = fgn;
    return
end
