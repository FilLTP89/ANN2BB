function [D1,nd1,Ia] = arias_intensity(x,dt,I1);
% x = acceleration time history 
% dt = sampling time step of x 
% I1 = target value of the arias intensity 
% output: 
% nd1 = position corresponding to I1 value of the arias intensity 
% D1 = time corresponding to I1 value of the arias intensity

gi = 9.81; %gravity acceleration
N=length(x);
Ia(1)=0;
for i=1:N-1,
    Ia(i+1)=Ia(i)+x(i+1)^2;
end
Ia = Ia.*pi/2/gi; 
Ia=Ia./max(Ia);
for i=2:N,
    if ((Ia(i)>=I1)&(Ia(i-1)<I1)),
        nd1=i;
    end
end
D1=(nd1)*dt;

return 

