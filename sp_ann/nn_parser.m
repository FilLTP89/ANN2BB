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
        ann.fn = varargin{2};
        load(ann.fn);
        ann.net = net;
    else
        % [TODO] generate ann
    end
    %%
    % _reference natural periods (last 2 values = PGV and PGD)_
    %----------------------------------------------------------------------
    TnC = varargin{1}; % corner period for ANN
    switch TnC
        case 0.5
            inp.vTn = [0.6:0.1:1.0,1.25:0.25:5.0]';
            tar.vTn = [0,0.05,0.1:0.1:0.5]';
        case  0.75
            inp.vTn = [0.8:0.1:1.0,1.25:0.25:5.0]';
            tar.vTn = [0,0.05,0.1:0.1:0.7,0.75]';
        case 1.0
            inp.vTn = (1.25:0.25:5.0)';
            tar.vTn = [0,0.05,0.1:0.1:1.0]';
    end
    inp.nT = numel(inp.vTn);
    tar.nT = numel(tar.vTn);
    ann.inp = inp;
    ann.tar = tar;
    ann.TnC = TnC;
    %% OUTPUT
    varargout{1} = ann;
    return
end