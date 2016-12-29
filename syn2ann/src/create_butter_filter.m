%% *Create Butterworth Filter*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _create_butterworth_filter_: function design a butterworth filter for
% acausal signal filtering in time (via filtfilt)
%% INPUT:
% _bfo (filter's order)_
% _lfr (corner frequency (for high-pass filtering))_
% _hfr (cut-off frequency (for low-pass filtering))_
% _fNy (Nyquist frequency)_
%% OUTPUT:
% _bfb (b coefficient of Butterworth's filter)
% _bfa (a coefficient of Butterworth's filter)
% _flag (true/false if Butterworth's filter was created or not)
function [varargout] = create_butter_filter(varargin)
    bfo = varargin{1};
    lfr = varargin{2};
    hfr = varargin{3};
    fNy = varargin{4};
    
    flag=true;
    
    if (~isempty(lfr))&&(~isempty(hfr))
        % BP filter definition
        [bfb,bfa] = butter(bfo,[lfr hfr]./fNy,'bandpass');
%         fprintf('\nbandpass filtering: %.2f-%.2f Hz\n',lfr,hfr);
    elseif (~isempty(lfr))&&(isempty(hfr))
        % LF filter definition
        [bfb,bfa] = butter(bfo,lfr./fNy,'high');
%          fprintf('\nhigh-pass filtering: %.2f Hz\n',lfr);
    elseif (isempty(lfr))&&(~isempty(hfr))
        % LF filter definition
        [bfb,bfa] = butter(bfo,hfr./fNy,'low');
%         fprintf('\nlow-pass filtering: %.2f Hz\n',hfr);
    elseif (isempty(lfr))&&(isempty(hfr))
        flag=false;
        bfb = -1;
        bfa = -1;
%         fprintf('\noriginal records\n');
    end
    
    varargout{1} = bfb;
    varargout{2} = bfa;
    varargout{3} = flag;
    return
end