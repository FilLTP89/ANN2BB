function [varargout] = spectral_scaling(varargin)
    %===============
    % Spectral Scaling
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % spectral_scaling: function to scale input response spectra so to
    % matcht the target spectra
    % INPUT:  inp_dtm (input signal time-step)
    %         inp_tha (input accelerogram)
    %         tar_psa (target psa response spectrum)
    %         tar_vTn (target natural periods)
    % OUTPUT: obj_dtm (matched record time-step)
    %         obj_tha (matched accelerogram)
    %         obj_thv (matched velocigram)
    %         obj_thd (matched displacement record)
    %         obj_frq (matched frequencies)
    %         obj_psa (matched psa response sepctrum)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    %----------------------------------------------------------------------
    % input/target spectra
    %----------------------------------------------------------------------
    inp_dtm = varargin{1};
    inp_tha = varargin{2};
    tar_psa = varargin{3};
    tar_vTn = varargin{4};
    tar_nTn = numel(tar_vTn);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % filtering properties
    %----------------------------------------------------------------------
    fac = 1;
    scl = 1;
    lfr = 0.2; % corner frequency
    hfr = 50;  % cutoff frequency
    frm = 40;  % max frequency considered
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % resampling and processing
    %----------------------------------------------------------------------
    [obj_dtm,obj_tha,obj_ntm,~] = seismo_rsmpl(inp_dtm,inp_tha,fac,scl);
    [obj_tha,~,~] = band_pass_filter(obj_tha,obj_dtm,lfr,hfr);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % spectral matching set up
    %----------------------------------------------------------------------
    nit = 5;                    % number of iterations
    obj_nfr = tar_nTn+2;
    obj_frq = -ones(obj_nfr,1);
    obj_frq(end)     = 0.5/obj_dtm;
    obj_frq(2:end-1) = flip(1./tar_vTn);
    %======================================================================
    %======================================================================
    % SPECTRAL MATCHING
    %======================================================================
    for i_ = 1:nit % spectral matching iterations
        %------------------------------------------------------------------
        % fourier spectra
        %------------------------------------------------------------------
        [~,obj_fsa,~,~,~,~] = super_fft(obj_dtm,obj_tha,0);
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % psa response spectra
        %------------------------------------------------------------------
        [~,~,~,obj_psa,~] = newmark_sd(obj_tha,obj_dtm,tar_vTn,0.05);
        obj_psa(1) = max(abs(obj_tha));
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % psa response spectral ratio
        %------------------------------------------------------------------
        obj_rra          = -ones(obj_nfr,1);
        obj_rra(1)       = 1.;
        obj_rra(end)     = obj_psa(1)./tar_psa(1);
        obj_rra(2:end-1) = flip(obj_psa./tar_psa);
        obj_rra = interp1(obj_frq,obj_rra,inp_frq,'linear');
        obj_rra(inp_frq>=frm) = obj_psa(1)./tar_psa(1);
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % new (scaled) fourier spectrum
        %------------------------------------------------------------------
        obj_fsa = obj_fsa./obj_rra;
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % new accelerogram
        %------------------------------------------------------------------
        obj_tha = detrend(real(ifft(obj_fsa),obj_ntm,'symmetric'),'constant');
        %------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
    % filtering
    %----------------------------------------------------------------------
    [obj_tha,obj_thv,obj_thd] = band_pass_filter(obj_tha,obj_dtm,lfr,hfr);
    %======================================================================
    varargout{1} = obj_dtm;
    varargout{2} = obj_tha;
    varargout{3} = obj_thv;
    varargout{4} = obj_thd;
    varargout{5} = obj_frq;
    varargout{6} = obj_psa;
return
end
