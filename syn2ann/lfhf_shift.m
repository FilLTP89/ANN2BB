%% *Shift time records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_shift_: shift hf records according to the difference in T5%
% between lf and hf records
%% INPUT:
% * _nss(low frequency records)_
% * _sps(high frequency records)_
%% OUTPUT:
% * _sps(high frequency records)_
function [varargout] = lfhf_shift(varargin)
    %% SET-UP
    nss = varargin{1};
    sps = varargin{2};
    %% ARIAS INTENSITY
    I1 = 0.05;
    [nss,sps] = lfhf_arias(nss,sps,I1);
    %% TIME SHIFTING HF
    for i_ = 1:nss.mon.na
        for j_ = 1:nss.mon.nc
            idx_lf = nss.syn{i_}.AI5.(nss.mon.cp{j_});
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
                    hf(k_) = sps.syn{i_}.tha.(sps.mon.cp{j_})(k_-abs(idx_df)+1);
                end
            end
            [sps.syn{i_}.tha.(nss.mon.cp{j_}),...
                sps.syn{i_}.thv.(nss.mon.cp{j_}),...
                sps.syn{i_}.thd.(nss.mon.cp{j_})] = ...
                idc_tha(sps.mon.dtm(i_),hf(:));
                
        end
    end
    %% OUTPUT
    varargout{1} = nss;
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
    nss = varargin{1};
    sps = varargin{2};
    I1 = varargin{3};
    %% ARIAS INTENSITY
    for i_ = 1:nss.mon.na
        for j_ = 1:nss.mon.nc
            % low frequency-arias intensity
            [~,nss.syn{i_}.AI5.(nss.mon.cp{j_}),...
                nss.syn{i_}.Ain.(nss.mon.cp{j_})] = ...
                arias_intensity(nss.syn{i_}.tha.(nss.mon.cp{j_}),nss.mon.dtm(i_),I1);
            % high frequency-arias intensity
            [~,sps.syn{i_}.AI5.(sps.mon.cp{j_}),...
                sps.syn{i_}.Ain.(sps.mon.cp{j_})] = ...
                arias_intensity(sps.syn{i_}.tha.(sps.mon.cp{j_}),sps.mon.dtm(i_),I1);
        end
    end
    %% OUTPUT
    varargout{1} = nss;
    varargout{2} = sps;
    return
end