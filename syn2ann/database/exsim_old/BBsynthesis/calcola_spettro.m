function [spettro,freq]=calcola_spettro(segnale,dt)
%funzione che calcola lo spettro di Fourier di un generico segnale:
%richiede come dati di input la variabile in cui è memorizzato il 
%segnale ed il suo passo temporale 

npun=length(segnale);

ig=0;
while (2^ig<=npun); 
    ig=ig+1;
end
Nfft=2^ig; %numero di punti su cui viene calcolata la FFT%
durata=(Nfft-1).*dt; %durata del segnale che include anche i punti fittizi introdotti per il calcolo della FFT%
df=1/durata;
freq=0:df:1/2/dt;

segnale(npun+1:Nfft)=0;
npun=Nfft;
t=[0:dt:(npun-1)*dt];
spettro=dt*fft(segnale,Nfft); %trasformata di Fourier dell'accelerogramma non corretto%
spettro = (spettro(1:length(freq))); 
return
