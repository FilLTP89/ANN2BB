% *Compute response spectra*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_pad_: function that pads the lf and hf records and compute
% arias intensity and T5% for lf and hf records to align hf to lf
%% INPUT:
% * _pbs(low frequency records)_
% * _sps(high frequency records)_
%% OUTPUT:
% * _pbs(low frequency records)_
% * _sps(high frequency records)_
function [varargout] = lfhf_pad(varargin)
    %% SET-UP
    pbs = varargin{1};
    sps = varargin{2};
    tpad = 0;
    for i_ = 1:pbs.mon.na
        %% PAD DEFINITION
        npad_lf = round(tpad/pbs.mon.dtm(i_));
        npad_hf = round(tpad/sps.mon.dtm(i_));
        npad_lf = pbs.mon.ntm(i_)+npad_lf;
        npad_hf = sps.mon.ntm(i_)+npad_hf;
        npad = max([npad_lf,npad_hf]);
        %% PADDING RECORDS
        for j_ = 1:pbs.mon.nc
            % tapering
            pbs.syn{i_}.tha.(pbs.mon.cp{j_}) = ...
                taper_fun(pbs.syn{i_}.tha.(pbs.mon.cp{j_}),5,0,1); 
            
            % _padding low-frequency_
            pbs.syn{i_}.tha.(pbs.mon.cp{j_}) = ...
                padarray(pbs.syn{i_}.tha.(pbs.mon.cp{j_}),...
                npad-pbs.mon.ntm(i_),0,'post');
            
            % _padding high-frequency_
            sps.syn{i_}.tha.(sps.mon.cp{j_}) = ...
                padarray(sps.syn{i_}.tha.(sps.mon.cp{j_}),...
                npad-sps.mon.ntm(i_),0,'post');
        end
        % _update time vector_
        pbs.mon.ntm(i_) = npad;
        sps.mon.ntm(i_) = npad;
        pbs.mon.vtm(i_) = {pbs.mon.dtm(i_)*(0:pbs.mon.ntm(i_)-1)'};
        sps.mon.vtm(i_) = {sps.mon.dtm(i_)*(0:sps.mon.ntm(i_)-1)'};
    end
    %% OUTPUT
    varargout{1} = pbs;
    varargout{2} = sps;
    return
end