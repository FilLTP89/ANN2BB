%% *Sabetta and Pugliese synthetics*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _arias intensity_: function to compute arias intensity and t95%
%% INPUT:
% * _tha (acceleration time history)_
% * _dtm (sampling time-step)_
% * _idx (target value of the arias intensity)_
%% OUTPUT:
% * _vtm_idx (time corresponding to idx value of the arias
% intensity)_
% * _idx (position corresponding to idx value of the arias
%         intensity)_
% * _Ain (Arias intensity)
function [varargout] = arias_intensity(varargin)
    %% *SET-UP*
    tha = varargin{1};
    dtm = varargin{2};
    idx = 0.05;
    if nargin>2
        idx = varargin{3};
    end
    ntm = length(tha);
    Ain = -ones(ntm,1);
    %
    % _gravity acceleration_
    %
    grv = 9.81;
    
    %% *ARIAS INTENSITY*
    Ain = cumtrapz(tha(2:end).^2);
    % Ain = cumsum(tha.^2);
    
    %
    % _normalization
    %
    Ain = (0.5.*pi.*dtm./grv).*Ain;
    Ain = Ain./max(Ain);
    %% *INDEX*
    idx = find(Ain<=idx,1,'last');
    vtm_idx = dtm*(idx-1);
    %% *OUTPUT*
    varargout{1} = vtm_idx;
    varargout{2} = idx;
    varargout{3} = Ain;
    return
end
