%% *Parse KKNPP record filename*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _parse_kknpp_file: function to parse and synchronize records from KKNPP database
%% INPUT:
% * fn (file name)
% * cp (motion directions)
%% OUTPUT:
% * dtm (time-step)
% * ntm (number of time-steps)
% * vtm (time-vector)
% * th (structure of parsed th)
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
       th.(cp{i_}) = data.(upper(cp{i_}))*scale;
    end
    ntm = numel(th.(cp{i_}));
    vtm = (0:ntm-1)*dtm;
    
    varargout{1} = dtm;
    varargout{2} = ntm;
    varargout{3} = vtm;
    varargout{4} = th;
    
    return
end