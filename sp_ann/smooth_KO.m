%================
%                    Konno Ohmachi Smooth function
%                       as = smooth_KO(af,b);
%
%  function  W = (sin(log10(f./fc)*b)./(log10(f./fc)).*b).^4;
%    ( fc: Central frequency [hz] around which the smoothing is performed)
%
% Input:        
%               af: function a(f)
%               b: bandwidth coefficient (usually b = 40)
%                  determines half-width of peak 10^(1/b)
%                   
% Output:       as: Smoothed function
%
% ! equally-spaced f
%-------------------------
% Silvana Montoya-Noguera (27/05/2013) 
%
%References:
% Konno, K. and Omachi, T., 1998, Bull. Seism. Soc. Am., 88, 228-241.
%===============

function as = smooth_KO(varargin)

af = varargin{1}; af = af(:);
b = varargin{2};

nt = size(af,1);
as(1:nt)=af;
fratio = 10^(2.5/b);

for ix = 2:nt
    ifc = ix-1; % center frequency index
    xl1 = max([ceil(ifc/fratio),1]);
    xl2 = min([fix(ifc*fratio),nt-1]);
    if xl2 >= nt-1
        xl2 = nt-1;
    end
    aas = b*log10((xl1:xl2)/(ix-1));
    aas(1:(ifc-xl1)) = (sin(aas(1:(ifc-xl1)))./aas(1:(ifc-xl1))).^4;
    aas((ifc-xl1+2):(xl2-xl1+1)) = (sin(aas((ifc-xl1+2):(xl2-xl1+1)))./aas((ifc-xl1+2):(xl2-xl1+1))).^4;
    aas(ifc-xl1+1) = 1;
    a1=sum((af((xl1+1):(xl2+1)))'.*aas);
    as(ix)=a1/sum(aas);
    clear a1
end;
