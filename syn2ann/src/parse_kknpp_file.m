%% *Parse KKNPP record filename*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _parse_kknpp_file_: function to parse and synchronize records from KKNPP database
%% *INPUT*
% * _fn (file name)_
%% *OUTPUT*
% * _dtm (time-step)_
% * _ntm (number of time-steps)_
% * _vtm (time-vector)_
% * _th (structure of parsed th)_
% * _mon_lon (station longitude)_
% * _mon_lat (station latitude)_

function [varargout] = parse_kknpp_file(varargin)
    %% *SET-UP*
    fn = varargin{1};
    cp = varargin{2};
    
    if ~ismatfile(fn)
        keyboard
        error('file extension not-valid');
    end
    scale = 1; % m/s2
    
    %% *READ DATA*
    data = load(fn);
    
    dtm = data.f^-1;
    for i_ = 1:numel(cp)
        if any(strcmpi(fieldnames(data),'Z'))
            if strcmpi(cp{i_}(1),'U')
                data.U=data.Z;
            end
        end
        th.(cp{i_}) = data.(upper(cp{i_}(1)))*scale;
    end
    ntm = numel(th.(cp{i_}));
    vtm = (0:ntm-1)*dtm;
    
    
    %% *OUTPUT*
    varargout{1} = dtm;
    varargout{2} = ntm;
    varargout{3} = vtm;
    varargout{4} = th;
    warning('check output!---> lat/lon added');
    varargout{5} = data.siteLon;
    varargout{6} = data.siteLat;
    varargout{7} = data.siteZ;
    varargout{8} = data.eventLon;
    varargout{9} = data.eventLat;
    
    return
end