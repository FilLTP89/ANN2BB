%% *Resampling records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_rsmpl_: function that resamples lf and hf records to have
% the same dt and number of points.
%% INPUT:
% * lf_sim (low frequency, numerical simulation records);
% * hf_sp96 (high frequency, Sabetta&Pugliese synthetic records);
% * dt (time step);
%% OUTPUT:
% * lf_sim (low frequency, numerical simulation records);
% * hf_sp96 (high frequency, Sabetta&Pugliese synthetic records);
% * lf_t (low frequency, numerical simulation time vector);
% * hf_t (high frequency, Sabetta&Pugliese synthetic time vector);
% * lf_dt (low frequency, numerical simulation time step);
% * hf_dt (high frequency, Sabetta&Pugliese synthetic time step);

function [varargout] = lfhf_rsmpl(varargin)
    %% SET-UP
    lf_sim = varargin{1};
    lf_t = varargin{2};
    hf_sp96 = varargin{3};
    hf_t = varargin{4};
    dt = varargin{5};
    lf_dt = dt;
    hf_dt = dt;
    
    %  _check time-step_
    [status,dtm,idx1,idx2] = check_dt(lf_dt,hf_dt);
    
    if status
    %% RESAMPLING
        for i_ = 1:numel(idx1)
            j_ = idx1(i_);
            k_ = idx2(i_);
            if logical(k_-1)
                fac = hf_dt(j_)/dtm(j_);
                [new_dtm,hf_sp96,new_ntm,new_vtm] = ...
                    structfun(@(v) seismo_rsmpl(hf_dt(j_),v,...
                    fac),hf_sp96,'UniformOutput',0);
                % length time vector (high frequency records)
                hf_nt(j_) = new_ntm;
                % time step (high frequency records))
                hf_dt(j_) = new_dtm;
                % time vector (high frequency records))
                hf_t(j_) = new_vtm;
            else
                fac = hf_dt(j_)/dtm(j_);
                [new_dtm,lf_sim,new_ntm,new_vtm] = ...
                    structfun(@(v) seismo_rsmpl(lf_dt(j_),v,...
                    fac),hf_sim,'UniformOutput',0);
                % length time vector (low frequency records)
                lf_nt(j_) = new_ntm;
                % time step (low frequency records))
                lf_dt(j_) = new_dtm;
                % time vector (low frequency records))
                lf_t(j_) = new_vtm;
            end
        end
    end
    
    if lf_dt == hf_dt
       dt = lf_dt;
    end
    varargout{1} = lf_sim;
    varargout{2} = lf_t;
    varargout{3} = hf_sp96;
    varargout{4} = hf_t;
    varargout{5} = dt;
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