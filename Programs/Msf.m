close all;
clear all;
pkg load signal;

#Read the data obtained from the sound card
f = fopen('../Data/GMDT_15s.dat', 'rb');
data = fread(f,inf,'int16'); %int16 pour short
fclose(f);
msf=data(2:3:end);
clear data;

msf=msf(3*192e3:end);
msf=msf-mean(msf);
msf=hilbert(msf);

##Characteristic variables
fs=192e3; %Sampling frequency
ts=1/fs; %Sampling period
freq=linspace(-fs/2,fs/2,length(msf)); %Frequency vector : [-fs/2:fs/2]
time=[0:ts:ts*(length(msf)-1)]; %Time vector

##Variables
tmp=0;

dec=120; %Decimating factor
fs_dec=fs/dec; %Sampling frequency before decimation
ts_dec=1/fs_dec; %Sampling period before decimation
freq_dec=linspace(-fs_dec/2,fs_dec/2,length(msf(1:120:end))); %Frequency vector before decimation
time_dec=[0:ts_dec:ts_dec*(length(msf(1:120:end))-1)]; %Time vector before decimation

%%%%%%%%%%%%%%%%%%%%%%MSF processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%____________________FIR 1_____________________%%
F=[0, 100, 200, fs/2]*2/fs; 

A=[1,  1,   0,   0];

b=firls(64,F,A);

%%____________________Transposition, filtering & decimating_____________________%%
fo=60e3;
lo=exp(j*2*pi*fo*time); %local oscillator

msf_shift=msf.*lo';

msf_shift=filter(b,1,msf_shift);
msf_shift=msf_shift(1:120:end);

##figure(1)
##  subplot(211)
##  plot(freq,fftshift(abs(fft(msf))));
##  title("MSF : Fourier transform of the signal acquired");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');
##  subplot(212)
##  plot(freq_dec,fftshift(abs(fft(msf_shift))));
##  title("MSF : Fourier transform of the signal acquired transposed and decimated");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');

%%____________________Transposition_____________________%%
f_b=find(freq_dec<=-500);f_b=f_b(end);
f_h=find(freq_dec>=500);f_h=f_h(1);

[~,fo]=max(fftshift(abs(fft(msf_shift))));
fo=freq_dec(fo);

lo=exp(j*2*pi*fo*time_dec); %Local oscillator
msf_shift=msf_shift.*lo';

%%____________________FIR 2_____________________%%
F=[0, 10, 30, fs_dec/2]*2/fs_dec; 

A=[1,  1,   0,   0];

b=firls(64,F,A);

msf_filt=filter(b,1,msf_shift);

##figure(2)
##subplot(211)
##plot(freq_dec,fftshift(abs(fft(msf_filt))));
##title("MSF : Fourier transform of the signal filtered");
##ylabel('Amplitude (a.u)');
##xlabel('Frequency (Hz)');
##subplot(212)
##plot(time_dec,abs(msf_filt));    %Abs to plot the real part of the signal
##title('MSF : amplitude evolution of the signal according to time');
##ylabel('Amplitude (a.u)');
##xlabel('Time (s)');

%%____________________Time datation_____________________%%
msf_t=abs(msf_filt);

threshold=0.6*mean(msf_t);
d_threshold=threshold*0.5/18;

figure(3)
plot(time_dec,msf_t);
##hold on;
##plot(time_dec,ones(1,length(msf_t))*threshold);
hold on,

msf_t2 = msf_t-threshold; 


Rank = find(msf_t2>=-d_threshold & msf_t2<=d_threshold);

##plot(time_dec(Rank), msf_t(Rank), 'x');

MRank = msf_t(Rank);   

%If points n & n+1 are less than 0.1s apart, we keep only the closest one to the threshold
for i=1:1:length(Rank)-1
  if (Rank(i+1)-Rank(i)<(0.1*1600))
    if (abs(threshold-(msf_t(Rank(i)+1)))<abs(threshold-(msf_t(Rank(i)))))
      Rank(i)=0;
    else
      Rank(i+1)=0;
    endif
  endif
endfor

tmp=1;
m_tmp=0;
for i=1:length(Rank)
  if (Rank(i)!=0)
    m_tmp(tmp)=Rank(i);
    tmp=tmp+1;
  endif
endfor
Rank=m_tmp;

hold on,

##plot(time_dec(Rank), msf_t(Rank), 'o');

%Selecting falling edge

for i=1:1:length(Rank)-1
  if (msf_t(Rank(i))<msf_t(Rank(i)+1))
    Rank(i)=0;
  endif
endfor

i=length(Rank);

if (msf_t(Rank(i))>msf_t(Rank(i)-1))
    Rank(i)=0;
endif

tmp=1;
m_tmp=0;
for i=1:length(Rank)
  if (Rank(i)!=0)
    m_tmp(tmp)=Rank(i);
    tmp=tmp+1;
  endif
endfor
Rank=m_tmp;


%We remove small errors

diff(Rank);


tmp=find(diff(Rank)<500);
Rank(tmp)=0;

for i=1:length(Rank)
  if (Rank(i)!=0)
    m_tmp(tmp)=Rank(i);
    tmp=tmp+1;
  endif
endfor
Rank=m_tmp;

hold on,
plot(time_dec(Rank), msf_t(Rank), 'd');

time_dec(Rank);

%%_______________________Time datation________________%%
lsp = 120*16; %1920 points

##figure(4)
for i=1:1:(length(Rank)-1)
  [u,v]=polyfit([Rank(i)-1:Rank(i)+1],msf_t([Rank(i)-1:Rank(i)+1])',1);
##  plot(v.yf); hold on;
  x_interp=linspace(Rank(i)-1,Rank(i)+1,lsp);
  y_interp=polyval(u,x_interp);
##  plot(y_interp),hold on,
  b=find(y_interp>=threshold);b=b(end);
  t_interp=linspace(time_dec(Rank(i)-1),time_dec(Rank(i)+1),lsp);
  rank(i)=(t_interp((b)-1));
endfor

index_secMSF=rank;