close all;
clear all;
pkg load signal;
graphique=1;

#Read the data obtained from the sound card
f = fopen('../Data/GMDT_15s.dat', 'rb');
data = fread(f,inf,'int16'); %int16 pour short
fclose(f);

gps=data(1:3:end);

gps=gps(5*192e3:end);

clear data;

##Characteristic variables
fs=192e3; %Sampling frequency
ts=1/fs; %Sampling period
freq=linspace(-fs/2,fs/2,length(gps)); %Frequency vector : [-fs/2:fs/2]
time=[0:ts:ts*(length(gps)-1)]; %Time vector

%%%%%%%%%%%%%%%%%%%%%%GPS processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
  plot(time, gps)
  xlabel('Time (s)');
  ylabel('Amplitude (a.u)');
  title('GPS : amplitude evolution of the signal according to time');

seuil_gps=0.5*max(gps);
k=find(gps<seuil_gps);gps(k)=0;
k=find(gps>seuil_gps);gps(k)=1;

tmp=find(diff(gps)==1);
index_secGPS=time(tmp);

