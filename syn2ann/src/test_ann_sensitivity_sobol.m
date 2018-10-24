%% *TEST ANN SENSITIVITY*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% Copyright 2018_
%% *NOTES*
% _test_ann_sensitivity_: function to test ANN sensitivity
%% *N.B.*
% Need for:

%% *REFERENCES*
% @Book{Book_Haykin_1999_ANN,
%   Title                    = {{Neural Newtworks: a Comprehensive Foundation}},
%   Author                   = {Haykin, S.},
%   Publisher                = {Prentice-Hall International, Inc.},
%   Year                     = {1999},
% 
%   File                     = {Book_Haykin_1999_ANN.pdf:Book_Haykin_1999_ANN.pdf:PDF}
% }

function [psa_out] = test_ann_sensitivity_sobol(varargin)
    trann_load_sensitivity;
    trann_test;
    return
end
