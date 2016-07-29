%% *Compute response/fourier spectra*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% sp_spectra: function to compute response/fourier spectra for a set
% of synthetics records generated with Sabetta and Pugliese method
%% INPUT:
%         _Tc (corner period)_
%         _fn (ann file name)_
%% OUTPUT:
%         _ann (artificial neural network structure)_
% N.B. Need for _newmark_sd.m_
%% REFERENCES:
% http://it.mathworks.com/help/nnet/ug/neural-network-object-properties.html?requestedDomain=www.mathworks.com

function [varargout] = nn_parser(varargin)
    
    %% SET-UP
    if nargin>1 % load existing ann
        ann = load(varargin{2},'net');
        ann.fn  = varargin{2};
    else
        % [TODO] generate ann
    end
    %%
    % _reference natural periods (last 2 values = PGV and PGD)_
    %----------------------------------------------------------------------
    TnC = varargin{1}; % corner period for ANN
    [inp.vTn,tar.vTn,inp.nT,tar.nT] = trann_define_inout(TnC);
    ann.inp = inp;
    ann.tar = tar;
    ann.TnC = TnC;
    %% OUTPUT
    varargout{1} = ann;
    return
end