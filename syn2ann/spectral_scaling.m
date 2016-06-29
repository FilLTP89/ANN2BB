%% *Spectral Scaling*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _spectral_scaling_: function to scale input response spectra so to
% matcht the target spectra
%% INPUT:
% * _inp_dtm (input signal time-step)_
% * _inp_tha (input accelerogram)_
% * _tar_psa (target psa response spectrum)_
% * _tar_vTn (target natural periods)_
%% OUTPUT:
% * _obj_dtm (matched record time-step)_
% * _obj_tha (matched accelerogram)_
% * _obj_thv (matched velocigram)_
% * _obj_thd (matched displacement record)_
% * _obj_vfr (matched frequencies)_
% * _obj_psa (matched psa response sepctrum)_
%% N.B.:
% Need for _band_pass_filter.m_
function [varargout] = spectral_scaling(varargin)
    
    %% SET-UP
    %%
    % _input/target spectra_
    inp_dtm = varargin{1};      % input time step
    inp_tha = varargin{2};      % input acceleration time history
    tar_psa = varargin{3};      % target psa
    tar_vTn = varargin{4};      % target natural periods
    tar_nT = numel(tar_vTn);   % number of target natural periods
    % _index of selected natual periods_
    if nargin>4
        tar_idx = varargin{5};
    else
        tar_idx = 1:tar_nT;
    end
    %%
    % _filtering properties_
    fac = 1;                    % resampling factor
    scl = 1;                    % scaling factor
    lfr = [];                  % corner frequency
    hfr = [];                   % cutoff frequency
    frm = 40;                   % max frequency considered
    %%
    % _spectral matching set up_
    obj_vfr = flip(1./tar_vTn); % target frequencies
    nit = 4;                    % number of iterations
    %%
    % _resampling_
    [obj_dtm,obj_tha,obj_ntm,~] = seismo_rsmpl(inp_dtm,inp_tha,fac,scl);
%     obj_dtm = inp_dtm;
%     obj_tha = inp_tha(:);
%     obj_ntm = numel(obj_tha);
    %%
    % _objective frequency_
    idx0 = find(tar_vTn == 0);  % indexes
    obj_vfr(tar_nT-idx0+1) = 0.5/obj_dtm;% Nyquist frequency
    tar_vTn(idx0) = 1e-3;
    obj_vfr = [0;obj_vfr];
    obj_vTn = tar_vTn;
    % _frequency vector_
    [inp_vfr,~,~,~,~,~] = ...
        super_fft(obj_dtm,obj_tha,0);
    idx = inp_vfr<0.5/obj_dtm;
    inp_vfr = inp_vfr(idx);
    obj_nfs = numel(inp_vfr);
    %% SPECTRAL MATCHING
    for i_ = 1:nit % spectral matching iterations
        %%
        % _fourier spectra_
        [~,~,~,obj_fsa,~,~] = super_fft(obj_dtm,obj_tha(:),0);
        obj_fsa=obj_fsa./obj_dtm;
        obj_fsa = obj_fsa(1:obj_nfs,1);
        obj_fsa(obj_nfs+1,1)  = 0.;
        obj_fsa(obj_nfs+2:2*obj_nfs,1) = flip(conj(obj_fsa(2:end-1,1)));
        %%
        % _psa response spectra at target periods_
        obj_psa = tar_psa(:);
        [~,~,~,obj_psa(tar_idx,1),~] = newmark_sd(obj_tha(:,1),obj_dtm,tar_vTn(tar_idx,1),0.05);
        % _constraint PGA_
        obj_psa(idx0,1) = max(abs(obj_tha));
        obj_psa = obj_psa(:,1);
        %%
        % _psa response spectral ratio_
        obj_rra = flip(obj_psa(:,1)./tar_psa(:,1));
        obj_rra = [1;obj_rra];
        obj_rra = interp1(obj_vfr(:,1),obj_rra(:,1),inp_vfr(:,1),'linear');
        %%
        % _constraint pga_
        obj_rra(inp_vfr>=frm,1) = obj_psa(1)./tar_psa(1);
        %%
        % _complex conjugate_
        obj_rra(obj_nfs+1,1)  = 1.;
        obj_rra(obj_nfs+2:2*obj_nfs,1) = flip(conj(obj_rra(2:end-1,1)));
        %%
        % _new (scaled) fourier spectrum_
        obj_fsa = obj_fsa./obj_rra;
        obj_fsa = obj_fsa(~isnan(obj_fsa));
        %%
        % _new accelerogram
        obj_tha = detrend(ifft(obj_fsa(:)));
        obj_tha = obj_tha(~isnan(obj_tha(1:obj_ntm)));
    end
    %%
    % _correct PGA on time history_
    % _filtering
    [obj_tha,obj_thv,obj_thd] = band_pass_filter(obj_dtm,obj_tha,lfr,hfr);
    obj_vfr = obj_vfr(2:end,1);
    
    %% OUTPUT
    varargout{1} = obj_dtm;
    varargout{2} = obj_tha;
    varargout{3} = obj_thv;
    varargout{4} = obj_thd;
    varargout{5} = obj_vfr;
    varargout{6} = obj_vTn;
    varargout{7} = obj_psa;
    return
end
