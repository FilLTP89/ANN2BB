%% *Shift time records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% align_lfhf : shift hf records according to the difference in T5%
% between lf and hf records
%% INPUT:
% * acc_sim, t_sim (low frequency records);
% * acc_sp96, t_sp96 (high frequency records);
%% OUTPUT:
% * acc_sp96,vel_sp96,dis_sp96(high frequency records);
function [varargout] = lfhf_align(varargin)
    %% SET-UP
    acc_sim = varargin{1};
    acc_sp96 = varargin{2};
    t = varargin{3};
    dt = varargin{4};
    t_sp96 = t;
    t_sim = t;
    
    %% ARIAS INTENSITY
    I1 = 0.05;
    [AI5_lf,AI5_hf] = lfhf_arias(acc_sim,acc_sp96,dt,I1);
    %% TIME SHIFTING HF
    
            idx_lf = AI5_lf;
            idx_hf = AI5_hf;
            idx_df = idx_hf-idx_lf;
            if idx_df>0
                cont = 1;
                for k_ = abs(idx_df):numel(t_sp96)
                    hf(cont) = acc_sp96(k_);
                    cont = cont +1;
                end
                for k_ = numel(t_sp96)+1:numel(t_sp96)+abs(idx_df)-1
                    hf(cont) = 0.0;
                    cont = cont +1;
                end
            else
                for k_ = 1:abs(idx_df)
                    hf(k_) = 0.0;
                end
                for k_ = abs(idx_df)+1:numel(t_sp96)
%                     hf(k_) = acc_sp96(k_-abs(idx_df)+1);
                    hf(k_) = acc_sp96(k_-abs(idx_df));
                end
            end
            
   acc_sp96 = hf(:);
   % velocity
   vel_sp96 = cumtrapz(acc_sp96)*dt;
   % displacement
   dis_sp96 = cumtrapz(vel_sp96)*dt;
                
    %% OUTPUT
    varargout{1} = acc_sp96;
    varargout{2} = vel_sp96;
    varargout{3} = dis_sp96;
    return
end
%% *Compute Arias Intensity*
% _Editor: Filippo Gatti/ modified by Maria Infantino
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _lfhf_arias_: compute arias intensity and T5% for lf and hf records
%% INPUT:
% * acc_sim, t_sim (low frequency records);
% * acc_sp96, t_sp96 (high frequency records);
% * I1 = 0.05
%% OUTPUT:
% * AI5_lf (position corresponding to 5% of the arias
%         intensity of low frequency records);
% * AI5_hf (position corresponding to 5% of the arias
%         intensity of high frequency records);
function [varargout] = lfhf_arias(varargin)
    %% SET-UP
    acc_sim = varargin{1};
    acc_sp96 = varargin{2};
    dt = varargin{3};
    I1 = varargin{4};
    %% ARIAS INTENSITY
      % low frequency-arias intensity
      [~,AI5_lf,Ain_lf] = arias_intensity(acc_sim,dt,I1);
      % high frequency-arias intensity
      [~,AI5_hf,Ain_hf] = arias_intensity(acc_sp96,dt,I1);
        
        
    %% OUTPUT
    varargout{1} = AI5_lf;
    varargout{2} = AI5_hf;
    return
end