close all;
clear all;
pkg load signal;
graphique=1;

#Read the data obtained from the sound card
f = fopen('../Data/GMDT_15s.dat', 'rb');
data = fread(f,inf,'int16'); %int16 pour short
fclose(f);
tdf=data(3:3:end);

clear data;

tdf=tdf(3*192e3:end);
tdf=tdf-mean(tdf);
tdf=hilbert(tdf);

##Characteristic variables
fs=192e3; %Sampling frequency
ts=1/fs; %Sampling period
freq=linspace(-fs/2,fs/2,length(tdf)); %Frequency vector : [-fs/2:fs/2]
time=[0:ts:ts*(length(tdf)-1)]; %Time vector

##Variables
tmp=0;
threshold = 0.6;
dec=120; %Decimating factor
fs_dec=fs/dec; %Sampling frequency before decimation
ts_dec=1/fs_dec; %Sampling period before decimation
freq_dec=linspace(-fs_dec/2,fs_dec/2,length(tdf(1:120:end))); %Frequency vector before decimation
time_dec=[0:ts_dec:ts_dec*(length(tdf(1:120:end))-1)]; %Time vector before decimation
Points_25ms=fs_dec/40; %Number of points equals to 25ms after decimation
f_fft=linspace(-fs/2,fs/2,1024); %Frequency vector for a FFT on 1024 points

%%%%%%%%%%%%%%%%%%%%%%TDF processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%____________________FIR 1_____________________%%
F=[0, 2000, 3000, fs/2]*2/fs; 

A=[1,  1,   0,   0];

b=firls(64,F,A);

%%____________________Transposition, filtering & decimating_____________________%%
fo=72e3;
lo=exp(j*2*pi*fo*time); %local oscillator

tdf_shift=tdf.*lo';

tdf_shift=filter(b,1,tdf_shift);
tdf_shift=tdf_shift(1:120:end);

##figure(1)
##  subplot(211)
##  plot(freq,fftshift(abs(fft(tdf))));
##  title("TDF : Fourier transform of the signal acquired");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');
##  subplot(212)
##  plot(freq_dec,fftshift(abs(fft(tdf_shift))));
##  title("TDF : Fourier transform of the signal acquired transposed and decimated");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');

%%____________________Transposition_____________________%%
f_b=find(freq_dec<=-500);f_b=f_b(end);
f_h=find(freq_dec>=500);f_h=f_h(1);

[~,fo]=max(fftshift(abs(fft(tdf_shift))));
fo=freq_dec(fo);

lo=exp(j*2*pi*fo*time_dec); %Local oscillator
tdf_shift=tdf_shift.*lo';
tdf_filt=tdf_shift;

##figure(2)
##  subplot(211)
##  plot(freq_dec,fftshift(abs(fft(tdf_filt))));
##  title("TDF : Fourier transform of the signal filtered");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');
##  subplot(212)
##  plot(time_dec,abs(tdf_filt));    %Abs to plot the real part of the signal
##  title('TDF : amplitude evolution of the signal according to time');
##  ylabel('Amplitude (a.u)');
##  xlabel('Time (s)');

%%_____________________Processing phase of TDF_____________________%%

tdf_p = angle(tdf_filt);
tdf_uw=unwrap(tdf_p)';

% [u,v]=polyfit(time_dec,tdf_uw,1);

##figure(3)
##  subplot(311)
##  plot(time_dec,tdf_p)
##  xlabel('Time (s)');
##  ylabel('Phase (rad)');
##  title("TDF : phase evolution of the signal according to time");
##  
##  subplot(312)
##  plot(time_dec,tdf_uw), hold on, plot(time_dec,v.yf)
##  xlabel('Time (s)');
##  ylabel('Phase (rad)');
##  title("TDF : phase evolution of the signal according to time (unwrapped)");
##  
##  subplot(313)
##  plot(time_dec,tdf_uw-v.yf)


i_ant=1;
for i=1:fs_dec:length(tdf_uw)
  tdf_detrend(i_ant:i)=detrend(tdf_uw(i_ant:i));
  i_ant=i+1;
endfor
tdf_detrend=[tdf_detrend,detrend(tdf_uw(i_ant:end))];

figure(3)
  subplot(311)
  plot(time_dec,tdf_p)
  xlabel('Time (s)');
  ylabel('Phase (rad)');
  title("TDF : phase evolution of the signal according to time");
  
  subplot(312)
  plot(time_dec,tdf_uw)
  xlabel('Time (s)');
  ylabel('Phase (rad)');
  title("TDF : phase evolution of the signal according to time (unwrapped)");
  
  subplot(313)
  plot(time_dec,tdf_detrend)
  xlabel('Time (s)');
  ylabel('Phase (rad)');
  title("TDF : phase evolution of the signal according to time (without frequency drift)");
  
%%________________Shifting, decimating & creating 3 possibles values________________%%

%Shifting on a high level
start=find(tdf_detrend(1:fs_dec) >= 0.8*max(tdf_detrend(1:2*fs_dec))); 
start=start(1);

tdf_shiftHL=tdf_detrend(start:end);

%Decimating the signal in such a way to have 1 point every 25ms
tmp=0;
tdf_40=tdf_shiftHL(1:Points_25ms:end);

%Creating 3 possibles values : -1 rad, 0rad et 1rad
k=find(tdf_40<-threshold);tdf_40(k)=-1;
k=find(tdf_40>threshold);tdf_40(k)=1;
k=find((tdf_40<threshold) & (tdf_40>-threshold));tdf_40(k)=0;

temps_2=[0:ts_dec*Points_25ms:ts_dec*Points_25ms*(length(tdf_40)-1)];

figure(4)
  title("Shifting and decimating")
  subplot(211)
  plot(time_dec(1:length(tdf_shiftHL)),tdf_shiftHL);
  title("TDF : phase evolution of the signal according to time (after shifting)");
  xlabel('Time (s)');
  ylabel('Phase (rad)');

  subplot(212)
  plot(temps_2,tdf_40)
  title("TDF : evolution de la phase au cours du temps (phase evolution of the signal according to time (after decimating and creating 3 possibles values)");
  xlabel('Time (s)');
  ylabel('Phase (rad)');

%%__________________Locating the pattern of the second____________%%
index_TDF=strfind(num2str(tdf_40'+1),num2str([0 0 0 0 1 0 -1 0]'+1));

tmp=1;
for k=2:length(index_TDF)
  if ((index_TDF(k) - index_TDF(k-1)) < 39) % *(k-tmp));
    index_TDF(k)=0;
    tmp=tmp+1;
  endif
endfor

tmp=1;
m_tmp=0;
for i=1:length(index_TDF)
  if (index_TDF(i)!=0)
    m_tmp(tmp)=index_TDF(i);
    tmp=tmp+1;
  endif
endfor
index_TDF=m_tmp;

%%_______________________Time datation________________%%
for i=1:(length(index_TDF)-1)
  o=(index_TDF(i)-1)*Points_25ms+start;
%  p=(index_TDF(i)+0)*Points_25ms+start;
%  q=(index_TDF(i)+1)*Points_25ms+start;
  a=(index_TDF(i)+2)*Points_25ms+start;
  b=(index_TDF(i)+3)*Points_25ms+start;
%  c=(index_TDF(i)+4)*Points_25ms+start;
  d=(index_TDF(i)+5)*Points_25ms+start;
  e=(index_TDF(i)+6)*Points_25ms+start;
  
  
  tdf_detrend(o:e)=tdf_detrend(o:e)-mean(tdf_detrend(o:a));
  
%  x_A=[0:1:b-a];
%  y_A=tdf_shiftHL(a:b);

%  x_B=[0:1:d-b];
%  y_B=tdf_shiftHL(b:d);

%  x_C=[0:1:e-d];
%  y_C=tdf_shiftHL(d:e);

%  [u_A,v_A]=polyfit(x_A,y_A,1);
%  [u_B,v_B]=polyfit(x_B,y_B,1);
%  [u_C,v_C]=polyfit(x_C,y_C,1);
  
  b_red=b+1;
  d_red=d-1;

  x_B_red=time_dec(b_red:d_red);   % /!\ temps
  y_B_red=tdf_detrend(b_red:d_red);
  [u_B_red,v_B_red]=polyfit(x_B_red,y_B_red,1);

%  code=[v_A.yf,v_B.yf,v_C.yf];

##   figure(5)
##   %Preuve du polyfit qui fonctionne + réduction
##   hold on, plot(tdf_detrend(o:e)), hold on,

   date=-u_B_red(2)/u_B_red(1);
   index_secTDF(i)=date; 
endfor