function [varargout] = nn_parser(varargin)
    %===============
    % Compute response/fourier spectra
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % sp_spectra: function to compute response/fourier spectra for a set
    % of synthetics records generated with Sabetta and Pugliese method
    % INPUT:  Tc (corner period)
    %         fn (ann file name)
    % OUTPUT: ann (artificial neural network structure)
    % N.B. Need for newmark_sd.m
    % REFERENCES: 
    % http://it.mathworks.com/help/nnet/ug/neural-network-object-properties.html?requestedDomain=www.mathworks.com
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    if nargin>1 % load existing ann
        ann.fn = varargin{2};
        load(ann.fn);
        ann.net = net;
    else
        % [TODO] generate ann
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % reference natural periods (last 2 values = PGV and PGD)
    %----------------------------------------------------------------------
    Tc = varargin{1}; % corner period for ANN
    switch Tc
        case 0.5
            inp.vTn = [0.6:0.1:1.0,1.25:0.25:5.0]';
            tar.vTn = [0,0.05,0.1:0.1:0.5]';
        case  0.75
            inp.vTn = [0.8:0.1:1.0,1.25:0.25:5.0]';
            tar.vTn = [0,0.05,0.1:0.1:0.7,0.75]';
        case 1.0
            inp.vTn = [1.25:0.25:5.0]';
            tar.vTn = [0,0.05,0.1:0.1:1.0]';
    end
    inp.nT = numel(inp.vTn);
    tar.nT = numel(tar.vTn);
    ann.inp = inp;
    ann.tar = tar;
    %==========================================================================
    varargout{1} = ann;
    return
end