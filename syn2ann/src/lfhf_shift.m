%% *Shift time records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_shift_: shift hf records according to the difference in T5%
% between lf and hf records
%% INPUT:
% * _pbs(low frequency records)_
% * _sps(high frequency records)_
%% OUTPUT:
% * _sps(high frequency records)_
function [varargout] = lfhf_shift(varargin)
    %% SET-UP
    pbs = varargin{1};
    sps = varargin{2};
    %% ARIAS INTENSITY
    I1 = 0.05;
    [pbs,sps] = lfhf_arias(pbs,sps,I1);
    %% TIME SHIFTING HF
    for i_ = 1:pbs.mon.na
        for j_ = 1:pbs.mon.nc
            idx_lf = pbs.syn{i_}.AI5.(pbs.mon.cp{j_});
            idx_hf = sps.syn{i_}.AI5.(sps.mon.cp{j_});
            idx_df = idx_hf-idx_lf;
            if idx_df>0
                cont = 1;
                for k_ = abs(idx_df):sps.mon.ntm(i_)
                    hf(cont) = sps.syn{i_}.tha.(sps.mon.cp{j_})(k_);
                    cont = cont +1;
                end
                for k_ = sps.mon.ntm(i_)+1:sps.mon.ntm(i_)+abs(idx_df)-1
                    hf(cont) = 0.0;
                    cont = cont +1;
                end
            else
                for k_ = 1:abs(idx_df)
                    hf(k_) = 0.0;
                end
                for k_ = abs(idx_df)+1:sps.mon.ntm(i_)
                    hf(k_) = sps.syn{i_}.tha.(sps.mon.cp{j_})(k_-abs(idx_df));
                end
            end
            [sps.syn{i_}.tha.(pbs.mon.cp{j_}),...
                sps.syn{i_}.thv.(pbs.mon.cp{j_}),...
                sps.syn{i_}.thd.(pbs.mon.cp{j_})] = ...
                idc_tha(sps.mon.dtm(i_),hf(:));
                
        end
    end
    %% OUTPUT
    varargout{1} = pbs;
    varargout{2} = sps;
    return
end
%% *Compute Arias Intensity*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_arias_: compute arias intensity and T5% for lf and hf records
%% INPUT:
% * _nss(low frequency records)_
% * _sps(high frequency records)_
% * I1 = 0.05
%% OUTPUT:
% * _nss(low frequency records)_
% * _sps(high frequency records)_
function [varargout] = lfhf_arias(varargin)
    %% SET-UP
    pbs = varargin{1};
    sps = varargin{2};
    I1 = varargin{3};
    %% ARIAS INTENSITY
    for i_ = 1:pbs.mon.na
        for j_ = 1:pbs.mon.nc
            % low frequency-arias intensity
            [~,pbs.syn{i_}.AI5.(pbs.mon.cp{j_}),...
                pbs.syn{i_}.Ain.(pbs.mon.cp{j_})] = ...
                arias_intensity(pbs.syn{i_}.tha.(pbs.mon.cp{j_}),pbs.mon.dtm(i_),I1);
            % high frequency-arias intensity
            [~,sps.syn{i_}.AI5.(sps.mon.cp{j_}),...
                sps.syn{i_}.Ain.(sps.mon.cp{j_})] = ...
                arias_intensity(sps.syn{i_}.tha.(sps.mon.cp{j_}),sps.mon.dtm(i_),I1);
        end
    end
    %% OUTPUT
    varargout{1} = pbs;
    varargout{2} = sps;
    return
end