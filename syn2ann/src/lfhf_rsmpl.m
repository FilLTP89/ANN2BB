%% *Resampling records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_rsmpl_: function that resamples lf and hf records to have
% the same dt and number of points.
%% INPUT:
% * _pbs (numerical simulation structure)_
% * _sps (Sabetta/Pugliese synthetics structure)_
%% OUTPUT:
% * _pbs (numerical simulation structure)_
% * _sps (Sabetta/Pugliese synthetics structure)_
function [varargout] = lfhf_rsmpl(varargin)
    %% SET-UP
    pbs = varargin{1};
    sps = varargin{2};
    %  _check time-step_
    [status,dtm,idx1,idx2] = check_dt(pbs.mon.dtm,sps.mon.dtm);
    
    if status
    %% RESAMPLING
        for i_ = 1:numel(idx1)
            j_ = idx1(i_);
            k_ = idx2(i_);
            if logical(k_-1)
                fac = sps.mon.dtm(j_)/dtm(j_);
                [new_dtm,sps.syn{j_}.tha,new_ntm,new_vtm] = ...
                    structfun(@(v) seismo_rsmpl(sps.mon.dtm(j_),v,...
                    fac),sps.syn{j_}.tha,'UniformOutput',0);
                sps.mon.ntm(j_) = new_ntm.(sps.mon.cp{1});
                sps.mon.dtm(j_) = new_dtm.(sps.mon.cp{1});
                sps.mon.vtm{j_} = new_vtm.(sps.mon.cp{1});
            else
                fac = pbs.mon.dtm(j_)/dtm(j_);
                [new_dtm,pbs.syn{j_}.tha,new_ntm,new_vtm] = ...
                    structfun(@(v) seismo_rsmpl(pbs.mon.dtm(j_),v,...
                    fac),pbs.syn{j_}.tha,'UniformOutput',0);
                
                pbs.mon.dtm(j_) = new_dtm.(sps.mon.cp{1});
                pbs.mon.ntm(j_) = new_ntm.(sps.mon.cp{1});
                pbs.mon.vtm{j_} = new_vtm.(sps.mon.cp{1});
            end
        end
    end
    
    varargout{1} = pbs;
    varargout{2} = sps;
    return
end
%% *Check time-steps and resample records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _check_dt_: function that check if time-step are similar and resamples
% the lf hf records accordingly.
%% INPUT:
% * _dt1 (time-step vector)_
% * _dt2 (time-step vector)_
%% OUTPUT:
% * _dtm (common time step)_
% * _ntm (common number of time-points)_
function [varargout] = check_dt(varargin)
    %% SET-UP
    status = 0;
    dt1 = varargin{1}(:);
    dt2 = varargin{2}(:);
    %% CHECK AND CORRECT DT
    [dtm,~ ] = min([dt1 dt2],[],2);
    [~,idx2] = max([dt1 dt2],[],2);
    idx1 = find(abs(dt1-dt2)>1e-10);
    
    if any(idx1)
        status  = 1;
    end
    %% OUTPUT
    varargout{1} = status;
    varargout{2} = dtm;
    varargout{3} = idx1;
    varargout{4} = idx2;
    return
end