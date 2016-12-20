% *Compute response spectra*
% _Editor: Filippo Gatti/modified by Maria Infantino
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_pad_: function that pads the lf and hf records;
%% INPUT:
% * lf_sim(low frequency records,SPEED simulations);
% * hf_sp96(high frequency records, sp96 synthetics);
% * dt(time step)
%% OUTPUT:
% * lf_sim(low frequency records,SPEED simulations);
% * hf_sp96(high frequency records, sp96 synthetic);
% * t(time vector)
function [varargout] = lfhf_pad(varargin)
    %% SET-UP
    lf_sim = varargin{1};
    hf_sp96 = varargin{2};
    dt = varargin{3};
    tpad = 0;
    
    %% PAD DEFINITION
    npad_lf = round(tpad/dt);
    npad_hf = round(tpad/dt);
    npad_lf = numel(lf_sim)+npad_lf;
    npad_hf = numel(hf_sp96)+npad_hf;
    npad = max([npad_lf,npad_hf]);
    %% PADDING RECORDS
    % tapering
    lf_sim = taper_fun(lf_sim,5,0,1); 
            
    % _padding low-frequency_
    lf_sim = padarray(lf_sim,npad-numel(lf_sim),0,'post');
            
    % _padding high-frequency_
    hf_sp96 = padarray(hf_sp96,npad-numel(hf_sp96),0,'post');
    
    % _update time vector_
    t = dt*(0:npad-1)';
    
    %% OUTPUT
    varargout{1} = lf_sim;
    varargout{2} = hf_sp96;
    varargout{3} = t;
    return
end