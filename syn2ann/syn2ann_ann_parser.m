%% *Compute response/fourier spectra*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_ann_parser_: function to parse trained ann networks
%% INPUT:
% * _Tc (corner period)_
% *_fn (ann file name)_
%% OUTPUT:
% *_ann (artificial neural network structure)_
% N.B. Need for _trann_define_inout.m_
%% REFERENCES:
% http://it.mathworks.com/help/nnet/ug/neural-network-object-properties.html?requestedDomain=www.mathworks.com

function [varargout] = syn2ann_ann_parser(varargin)
    %% *SET-UP*
    ann     = load(varargin{2},'net');
    ann.fn  = varargin{2};
    ann.TnC = varargin{1};
    
    %% *DEFINE INPUT/TARGET PERIODS*
    [ann.inp.vTn,ann.tar.vTn,ann.inp.nT,ann.tar.nT] = ...
        trann_define_inout(ann.TnC);
    
    %% OUTPUT
    varargout{1} = ann;
    return
end
