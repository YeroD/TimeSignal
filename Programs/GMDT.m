close all;
clear all;
pkg load signal;
graphique=1;

#Read the data obtained from the sound card
f = fopen('../Data/GMDT_15s.dat', 'rb');
data = fread(f,inf,'int16'); %int16 pour short
fclose(f);

gps=data(1:3:end);
msf=data(2:3:end);
tdf=data(3:3:end);

gps=gps(5*192e3:end);
msf=msf(5*192e3:end);
tdf=tdf(5*192e3:end);

clear data;

msf=msf-mean(msf);
msf=hilbert(msf);
dcf=msf;

tdf=tdf-mean(tdf);
tdf=hilbert(tdf);

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
Points_25ms=fs_dec/40; %Number of points equals to 25ms after decimation

%%%%%%%%%%%%%%%%%%%%%%GPS processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
  subplot(511)
  plot(time, gps)
  xlabel('Time (s)');
  ylabel('Amplitude (a.u)');
  title('GPS : amplitude evolution of the signal according to time');

seuil_gps=0.5*max(gps);
k=find(gps<seuil_gps);gps(k)=0;
k=find(gps>seuil_gps);gps(k)=1;

tmp=find(diff(gps)==1);
index_secGPS=time(tmp);

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

##figure(2)
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

##figure(3)
##  subplot(211)
##  plot(freq_dec,fftshift(abs(fft(msf_filt))));
##  title("MSF : Fourier transform of the signal filtered");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');
##  subplot(212)
##  plot(time_dec,abs(msf_filt));    %Abs to plot the real part of the signal
##  title('MSF : amplitude evolution of the signal according to time');
##  ylabel('Amplitude (a.u)');
##  xlabel('Time (s)');

figure(1)
  subplot(512)
  plot(time_dec,abs(msf_filt));    %Abs to plot the real part of the signal
  title('MSF : amplitude evolution of the signal according to time');
  ylabel('Amplitude (a.u)');
  xlabel('Time (s)');;

%%____________________Locating the pattern of the second_____________________%%
msf_t=abs(msf_filt);

threshold=0.6*mean(msf_t);
d_threshold=threshold*0.5/18;

##figure(4)
##plot(time_dec,msf_t);
##hold on;
##plot(time_dec,ones(1,length(msf_t))*threshold);
##hold on,

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

##hold on,

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

##figure(5)
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

%%%%%%%%%%%%%%%%%%%%%%DCF processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%____________________FIR_____________________%%
F=[0, 1000, 1500, fs/2]*2/fs; 

A=[1,  1,   0,   0];

b=firls(64,F,A);

%%____________________Transposition, filtering & decimating_____________________%%
fo=77.5e3;
lo=exp(j*2*pi*fo*time); %local oscillator

dcf_shift=dcf.*lo';

dcf_shift=filter(b,1,dcf_shift);
dcf_shift=dcf_shift(1:120:end);

##figure(6)
##  subplot(211)
##  plot(freq,fftshift(abs(fft(dcf))));
##  title("DCF : Fourier transform of the signal acquired");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');
##  subplot(212)
##  plot(freq_dec,fftshift(abs(fft(dcf_shift))));
##  title("DCF : Fourier transform of the signal acquired transposed and decimated");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');
##return
%%____________________Transposition_____________________%%
[~,fo]=max(fftshift(abs(fft(dcf_shift))));
fo=freq_dec(fo);

lo=exp(j*2*pi*fo*time_dec); %Local oscillator
dcf_shift=dcf_shift.*lo';

dcf_filt=dcf_shift;

##figure(7)
##  subplot(211)
##  plot(freq_dec,fftshift(abs(fft(dcf_filt))));
##  title("DCF : Fourier transform of the signal filtered");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');
##  subplot(212)
##  plot(time_dec,abs(dcf_filt));    %Abs to plot the real part of the signal
##  title('DCF : amplitude evolution of the signal according to time');
##  ylabel('Amplitude (a.u)');
##  xlabel('Time (s)');

%%____________________FIR 2_____________________%%
F=[0, 50, 100, fs_dec/2]*2/fs_dec; 

A=[1,  1,   0,   0];

b=firls(64,F,A);

figure(1)
  subplot(514)
  plot(time_dec,abs(filter(b,1,dcf_filt))); 
##  plot(time_dec,abs(dcf_filt)); 
  xlabel('Time (s)');
  ylabel('Amplitude (a.u)');
  title('DCF : amplitude evolution of the signal according to time');
  
%%____________________Locating the pattern of the second_____________________%%
dcf_t=abs(filter(b,1,dcf_filt));

threshold=0.6*mean(dcf_t);
d_threshold=threshold/9;

##figure(4)
##plot(time_dec,dcf_t);
##hold on;
##plot(time_dec,ones(1,length(dcf_t))*threshold);
##hold on,

dcf_t2 = dcf_t-threshold; 


Rank = find(dcf_t2>=-d_threshold & dcf_t2<=d_threshold);

##plot(time_dec(Rank), dcf_t(Rank), 'x');

MRank = dcf_t(Rank);   

%If points n & n+1 are less than 0.1s apart, we keep only the closest one to the threshold
for i=1:1:length(Rank)-1
  if (Rank(i+1)-Rank(i)<(0.1*1600))
    if (abs(threshold-(dcf_t(Rank(i)+1)))<abs(threshold-(dcf_t(Rank(i)))))
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

##plot(time_dec(Rank), dcf_t(Rank), 'o');

%Selecting falling edge

for i=1:1:length(Rank)-1
  if (dcf_t(Rank(i))<dcf_t(Rank(i)+1))
    Rank(i)=0;
  endif
endfor

i=length(Rank);

if (dcf_t(Rank(i))>dcf_t(Rank(i)-1))
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
  plot(time_dec(Rank), dcf_t(Rank), 'd');

time_dec(Rank);

%%____________________Phase demodulation_____________________%%
dcf_p = angle(dcf_filt);

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

##figure(8)
##  subplot(311)
##  plot(time_dec,dcf_p)
##  title('DCF : phase evolution of the signal according to time');
##  xlabel('Time (s)');
##  ylabel('Phase (rad)');
##
##  subplot(312)
##  plot(time_dec,dcf_uw)
##  title('DCF : phase evolution of the signal according to time (withouth phase drift)');
##  xlabel('Time (s)');
##  ylabel('Phase (rad)');
##
##  subplot(313)
##  plot(time_dec,abs(correlation)), 
##  title('DCF : correlation between the prn code and phase');
##  xlabel('Time (s)');
##  ylabel('Amplitude (a.u)');
  
figure(1)
  subplot(515)
  plot(time_dec,abs(correlation)), 
  xlabel('Time (s)');
  ylabel('Amplitude (a.u)');
  title('DCF : correlation between the prn code and phase');

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

##figure(9)
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

##figure(10)
##  plot(freq_dec,fftshift(abs(fft(tdf_filt))));
##  title("TDF : Fourier transform of the signal filtered");
##  ylabel('Amplitude (a.u)');
##  xlabel('Frequency (Hz)');

%%_____________________Processing phase of TDF_____________________%%

tdf_p = angle(tdf_filt);
tdf_uw=unwrap(tdf_p)';

% [u,v]=polyfit(time_dec,tdf_uw,1);

##figure(11)
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

##figure(12)
##  subplot(311)
##  plot(time_dec,tdf_p)
##  xlabel('Time (s)');
##  ylabel('Phase (rad)');
##  title("TDF : phase evolution of the signal according to time");
##  
##  subplot(312)
##  plot(time_dec,tdf_uw)
##  xlabel('Time (s)');
##  ylabel('Phase (rad)');
##  title("TDF : phase evolution of the signal according to time (unwrapped)");
##  
##  subplot(313)
##  plot(time_dec,tdf_detrend)
##  xlabel('Time (s)');
##  ylabel('Phase (rad)');
##  title("TDF : phase evolution of the signal according to time (without frequency drift)");
  
%%________________Shifting, decimating & creating 3 possibles values________________%%

%Shifting on a high level
start=find(tdf_detrend(1:fs_dec) >= 0.8*max(tdf_detrend(1:2*fs_dec))); 
start=start(1);

tdf_shiftHL=tdf_detrend(start:end);

%Decimating the signal in such a way to have 1 point every 25ms
tmp=0;
tdf_40=tdf_shiftHL(1:Points_25ms:end);

threshold=0.6;
%Creating 3 possibles values : -1 rad, 0rad et 1rad
k=find(tdf_40<-threshold);tdf_40(k)=-1;
k=find(tdf_40>threshold);tdf_40(k)=1;
k=find((tdf_40<threshold) & (tdf_40>-threshold));tdf_40(k)=0;

temps_2=[0:ts_dec*Points_25ms:ts_dec*Points_25ms*(length(tdf_40)-1)];

##figure(13)
##  title("Shifting and decimating")
##  subplot(211)
##  plot(time_dec(1:length(tdf_shiftHL)),tdf_shiftHL);
##  title("TDF : phase evolution of the signal according to time (after shifting)");
##  xlabel('Time (s)');
##  ylabel('Phase (rad)');
##
##  subplot(212)
##  plot(temps_2,tdf_40)
##  title("TDF : phase evolution of the signal according to time (after decimating and creating 3 possibles values)");
##  xlabel('Time (s)');
##  ylabel('Phase (rad)');

figure(1)
  subplot(513)
  plot(time_dec,tdf_detrend)
  xlabel('Time (s)');
  ylabel('Phase (rad)');
  title("TDF : phase evolution of the signal according to time");
  
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

  hold on,
  pos=(index_TDF+4)*Points_25ms+start;
  plot(time_dec(pos), tdf_detrend(pos), 'd');

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

##   figure(14)
##   %Preuve du polyfit qui fonctionne + réduction
##   hold on, plot(tdf_detrend(o:e)), hold on,

   date=-u_B_red(2)/u_B_red(1);
   index_secTDF(i)=date; 
endfor