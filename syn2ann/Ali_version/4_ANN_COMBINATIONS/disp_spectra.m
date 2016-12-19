function ymax=disp_spectra(u,dt,Tn,zeta)
% calcola spettro di spostamento per un accelerogramma (u), letto come
% una colonna di punti equispaziati con passo dt
npun=length(u);
if (Tn==0),
    On=1e9;
else
    On=2*pi/Tn;
end;

fac1=1/dt^2;
fac2=zeta*On/dt;

%condizioni iniziali di quiete
ym1=0;
y0=0;
y=ones(1,npun);
y(1)=ym1;
y(2)=y0;

for i=2:npun-1,
   yp1=((2*fac1-On^2)*y0+(-fac1+fac2)*ym1-u(i))/(fac1+fac2);
   y(i+1)=yp1;
   ym1=y0;
   y0=yp1;
end

ymax=max(abs(y));
