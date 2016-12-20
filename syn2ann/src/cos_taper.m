% cos_taper.m
% Applies a 10% cosine taper
% usage:
% tapered = cos_taper(data);
function [varargout] = cos_taper(data)
    data = data(:);
    Nd   = length(data);
%     M=floor(Nd/40+0.5);
    M = round(Nd/40);
    
    idx0 = (M+1:Nd-M)';
    idx1 = (0:M-1)';
    idx2 = (Nd-M+1:Nd)';
    
    tapered         = zeros(Nd,1);
    tapered(idx0)   = data(idx0);
    tapered(idx1+1) = data(idx1+1).*(0.5.*(1+cos(flip(idx1).*pi./(M-1))));
    tapered(idx2)   = data(idx2).*(0.5.*(1+cos(idx1.*pi./(M-1))));
%     tw = tukeywin(Nd,0.05);
%     tapered = data(:).*tw(:);
    varargout{1} = tapered(:);
    return
end
