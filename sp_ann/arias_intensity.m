function [varargout] = arias_intensity(varargin)
    % INPUT:  tha (acceleration time history)
    %         dtm (sampling time-step)
    %         idx (target value of the arias intensity)
    % OUTPUT: vtm_idx (time corresponding to idx value of the arias intensity)
    %         idx (position corresponding to idx value of the arias intensity)
    %         Ain (Arias intensity)
    %======================================================================
    % SET-UP
    %======================================================================
    tha = varargin{1};
    dtm = varargin{2};
    idx = varargin{3};
    ntm = length(tha);
    Ain = -ones(ntm,1);
    %----------------------------------------------------------------------
    % gravity acceleration
    %----------------------------------------------------------------------
    grv = 9.81;
    %======================================================================
    %======================================================================
    % ARIAS INTENSITY
    %======================================================================
    Ain(1)=0;
%     Ain(2:end) = cumsum(tha(2:end).^2);
    Ain(2:end) = cumtrapz(tha(2:end).^2);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % normalization
    %----------------------------------------------------------------------
    Ain = (0.5.*pi.*dtm./grv).*Ain;
    Ain = Ain./max(Ain);
    %======================================================================
    %======================================================================
    % INDEX
    %======================================================================
    idx = find(Ain>=idx,1,'first');
    vtm_idx = dtm*idx;
    %======================================================================
    varargout{1} = vtm_idx;
    varargout{2} = idx;
    varargout{3} = Ain;
    return
end
