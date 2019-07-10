close all;
clear all;
pkg load signal;
graphique=1;

#Read the data obtained from the sound card
f = fopen('../Data/GMDT_15s.dat', 'rb');
data = fread(f,inf,'int16'); %int16 pour short
fclose(f);
dcf=data(2:3:end);
clear data;

dcf=dcf(3*192e3:end);
dcf=dcf-mean(dcf);
dcf=hilbert(dcf);

##Characteristic variables
fs=192e3; %Sampling frequency
ts=1/fs; %Sampling period
freq=linspace(-fs/2,fs/2,length(dcf)); %Frequency vector : [-fs/2:fs/2]
time=[0:ts:ts*(length(dcf)-1)]; %Time vector
N=128; %Convolution windows


##Variables
tmp=0;
threshold = 0.6;
dec=120; %Decimating factor
fs_dec=fs/dec; %Sampling frequency before decimation
ts_dec=1/fs_dec; %Sampling period before decimation
freq_dec=linspace(-fs_dec/2,fs_dec/2,length(dcf)); %Frequency vector before decimation
time_dec=[0:ts_dec:ts_dec*(length(dcf(1:120:end))-1)]; %Time vector before decimation
Points_25ms=fs_dec/40; %Number of points equals to 25ms after decimation
f_fft=linspace(-fs/2,fs/2,1024); %Frequency vector for a FFT on 1024 points

%%%%%%%%%%%%%%%%%%%%%%DCF processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%____________________Variable_____________________%%
f_b=find(f_fft<=77000);f_b=f_b(end);
f_h=find(f_fft>=78000);f_h=f_h(1);

[~,fo]=max(fftshift(abs(fft(dcf,1024)))(f_b:f_h));
fo=f_fft(fo+f_b-1);

%%____________________Filtering_____________________%%
lo=exp(j*2*pi*fo*time); %Local oscillator

dcf_shift=dcf.*lo';
filtre_dcf=conv(dcf_shift,ones(N,1)/N); %Filter to keep only DCF frequency
filtredcf_fin=filtre_dcf(N/2:end-N/2);

##figure(1)
##plot(freq,fftshift(abs(fft(dcf))));
##title("DCF : Fourier transform of the signal acquired");
##hold on;
##plot(freq,fftshift(abs(fft(dcf_shift))), 'r');
##hold on;
##plot(freq,fftshift(abs(fft(filtredcf_fin))), 'g');
##ylabel('Amplitude (a.u)');
##xlabel('Frequency (Hz)');
##legend('Signal acquired','Signal transposed');
 

filtredcf_fin=filtredcf_fin(1:120:end);
f_tmp=linspace(-fs/2,fs/2,length(filtredcf_fin)); %Frequency vector
t_tmp=[0:ts:ts*(length(filtredcf_fin)-1)]; %Time vector

%plot(f_tmp,fftshift(abs(fft(filtredcf_fin)))),

[~,fo]=max(fftshift(abs(fft(filtredcf_fin))));
fo=f_tmp(fo);

lo=exp(j*2*pi*fo*t_tmp); %local oscillator
filtredcf_fin=filtredcf_fin.*lo';

dcf_t=abs(filtredcf_fin);
k=find(dcf_t<19);dcf_t(k)=0;
k=find(dcf_t>=19);dcf_t(k)=1;



##figure(1)
##plot(time_dec,dcf_t)



%plot(f_tmp,fftshift(abs(fft(filtredcf_fin)))), return,

%%____________________Phase demodulation_____________________%%
dcf_p = angle(filtredcf_fin);

i_ant=1;
for i=1:fs:length(dcf_p)
  dcf_alt(i_ant:i)=detrend(dcf_p(i_ant:i));
  i_ant=i+1;
endfor
dcf_alt=[dcf_alt,detrend(dcf_p(i_ant:end))'];
dcf_uw=dcf_alt;


%%____________________Correlation_____________________%%
A=floor(fs_dec/(77500/120));

prn = load('prn.dat');

increment=192000*120/77500/120; 
courant=0.1;
prni=[];
for p=1:512
  prni=[prni ones(1,(floor(courant+increment)-ceil(courant))+1)*prn(p)];
  courant=courant+increment;
endfor
prn=prni;


prn = prn - mean(prn);

correlation = xcorr(dcf_p,prn);
correlation = correlation(floor(length(correlation)/2)+1:end);

abs_corr=abs(correlation);

figure(2)
 subplot(311)
 plot(time_dec,dcf_p)
 title('DCF : phase evolution of the signal according to time');
 xlabel('Time (s)');
 ylabel('Phase (rad)');

 subplot(312)
 plot(time_dec,dcf_uw)
 title('DCF : phase evolution of the signal according to time (withouth phase drift)');
 xlabel('Time (s)');
 ylabel('Phase (rad)');

 subplot(313)
 plot(time_dec,abs(correlation)), 
 title('DCF : correlation between the prn code and phase');
 xlabel('Time (s)');
 ylabel('Amplitude (a.u)');

%%_______________________Time datation________________%%
[~,pic]=max(abs_corr(1:fs_dec));

place=1;
for i=pic:fs_dec:length(abs_corr)
  [~,rank]=max(abs_corr(i-Points_25ms:i+Points_25ms));
  m_tmp(place)=time_dec(rank+i-Points_25ms-1);
  Rank(place)=rank+i-Points_25ms-1;
  place=place+1;
endfor

m_tmp=0; tmp=1;
lsp=120*16; % 1920 points

for i=1:1:length(Rank)
  [u,v]=polyfit([Rank(i)-1:Rank(i)+1],abs_corr([Rank(i)-1:Rank(i)+1])',2);
##  plot(v.yf); hold on;
  x_interp=linspace(Rank(i)-1,Rank(i)+1,lsp);
  y_interp=polyval(u,x_interp);
##  plot(y_interp),hold on,
  [a,b]=max(y_interp);
  t_corr=linspace(time_dec(Rank(i)-1),time_dec(Rank(i)+1),lsp);
  rank(i)=(t_corr((b)-1));
endfor

index_secDCF=sort(rank);