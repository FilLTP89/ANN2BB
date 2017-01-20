function diff1_a=freq_differentiate_new(a,sps)
    dt = 1/sps;
% close all
% clc
% clear all
% 
% % a generic random waveform to be differentiated
% dt=0.01;f1=0.25;
% w1=2*pi*f1;t=0:dt:10;
% TW=tukeywin(length(t),0.5);
% for i=1:10
%     multA(i)=random('uni',-1,1);
% end
% a=zeros(1,length(TW));
% for i=1:length(t)
%     for j=1:10
%         a(i)=a(i)+multA(j)*sin(w1*j*t(i));
%     end
%     a(i)=a(i)*TW(i);
% end
% 
% a2=[zeros(1,length(a)) a zeros(1,length(a))];
% 
% clear a, clear t
% 
% a=a2;
% 
% clear a2
% 
t=[0:1:length(a)-1]*dt;max_t=max(t);
% %

% frequency content
N=2^nextpow2(length(t));
A=fft(a,N);
df=1/(N*dt);

f=[0:1:(N-1)]*df;
%

% frequency domain differentation
for i=1:N
    if i<=N/2+1
        f1(i)=f(i);
    else
        f1(i)=-f(N-i+2);
    end
    diff_A(i)=2*pi*f1(i)*sqrt(-1)*A(i);
    t1(i)=(i-1)*dt;
end
        
diff1_a=real(ifft(diff_A(:)));
% %
% 
% 
% %classical central difference
% diff2_a=zeros(1,length(a));
% 
% for i=2:length(t)-1
%     diff2_a(i)=(a(i+1)-a(i-1))/(2*dt);
% end
% %
% 
% 
% figure
% plot(t,a)
% 
% figure
% plot(t1,diff1_a,'linewidth',2,'color','black'),hold on
% plot(t,diff2_a,'linewidth',1,'color','red'),hold on
% xlim([0 max_t]);
% 
% 
% 
