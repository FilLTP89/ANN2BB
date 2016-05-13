function [varargout] = vel2acc(varargin)
    % ===============
    % Signal Processing
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % Politecnico di Milano - DICA
    % Copyright 2016
    % NOTES
    % vel2acc: function to differentiate velocigram to get accelerogram 
    % by applying butterworth filter and detrending.
    % INPUT:  dt (sampling time step)
    %         vel (input velocigram)
    %         fcut_low (corner frequency)
    %         fcut_high (cut-off frequency)
    % OUTPUT: acc (band-pass filtered acceleration time-history column vector)
    %         vel (velocity time-history column vector)
    %         dis (displacement time-history column vector)
    % ===============
    dt=varargin{1};                     % time step
    vel=varargin{2}(:);                 % velocigram
    fcut_low  =.01;                     % default corner frequency
    fcut_high = 25;                     % default cutoff frequency
    if nargin>=3
        fcut_low=varargin{3};           % customized corner frequency
    end
    if nargin>=4
        fcut_high=varargin{4};          % customized corner frequency
    end
    %======================================================================
    % BUTTERWORTH FILTER DEFINITION
    %======================================================================
    if nargin>3 || nargin<3
        % BP filter definition
        [bb,ab]=butter(2,[fcut_low fcut_high]*2*dt,'bandpass');
        fprintf('\nBP FILTER: f(corner): %.2f Hz - f(cut-off): %.2f Hz\n',...
            fcut_low,fcut_high);
    else
        % LF filter definition
        [bb,ab]=butter(2,fcut_low*2*dt,'high');
        fprintf('\nLF FILTER: f(corner): %.2f Hz\n',fcut_low);
    end
    %======================================================================
    % INTEGRATION/DIFFERENTIATION/DETRENDING/FILTERING
    %======================================================================
    vel=filtfilt(bb,ab,vel);
    % detrending velocity...
    vel=detrend(vel);
    % integrating displacement...
    dis=cumtrapz(vel)*dt;
    % filtering displacement...
    dis=filtfilt(bb,ab,dis);
    % detrending displacement...
    dis=detrend(dis);
    % cosinus-tapering displacement...
    dis=cos_taper(dis);
    % differentiating...
    vel=[0;diff(dis)/dt];
    acc=[0;diff(vel)/dt];
    
    varargout{1} = acc;
    varargout{2} = vel;
    varargout{3} = dis;
    return
end