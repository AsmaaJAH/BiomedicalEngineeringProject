close all;
clear all;
clc;
load("ecg.mat");
fs=500;
%sampling
t=1:1:length(ecg);
T=0:1/fs:length(ecg);
figure(1)
plot(T(t),ecg(t));
title('ecg signal');

%******0.5HZ*****%
% fourier transform
s(t)=fft(ecg(t));
s(t)=fftshift(s(t));
%filter****%
filter1=ones(1,(4170/2)-4);

filter3=zeros(1,8);
filter=[filter1 filter3 filter1] ;
% figure;
% plot(filter);
g(t)=s(t).*filter(t);
%inverse fourier transform
ecg1(t)=real(ifft(ifftshift(g(t))));

figure(2)
plot(T(t),ecg1(t),'LineWidth',0.1);
title('filter at 0.5HZ');
%********notch 50HZ*********/
%fourier transform

filter1=ones(1,(4170/2)+413);
filter2=ones(1,(4170/2)-421);
filter3=zeros(1,8);% figure;
% plot(filter);
filter=[filter1 filter3 filter2] ;
% figure;
% plot(filter);
ecg2(t)=ecg1(t).*filter(t);
%inverse fourier transform

figure(3)
plot(T(t),ecg2(t),'LineWidth',0.1);
title('notch filter at 50HZ');
%******************%
figure(4)
% fourier transform
fc_vec_sample=[41 , 166 , 333,500];
fc_HZ=[5 ,20 ,40,60];

for i=1:1:4
s3(t)=fft(ecg2(t));
s3(t)=fftshift(s3(t));
filter1=zeros(1,(4170/2)-fc_vec_sample(i));
filter3=ones(1,2*fc_vec_sample(i));
%as compromise on the cutt f freq we choose 50 HZ filter
%to not details in the ecg without give more noise
filter=[filter1 filter3 filter1] ;
% figure;
% plot(filter);
g3(t)=s3(t).*filter(t);

%inverse fourier transform
ecg3_(t)=real(ifft(ifftshift(g3(t)),'symmetric'));

subplot(4,1,i);
plot(T(t),ecg3_(t),'LineWidth',0.1);
title( [num2str(fc_HZ(i)),'HZ']);
end
%**************************** signal after increasing SNR**********
%WE CHOOSE 50 HZ as a cut off frequency
s3(t)=fft(ecg2(t));
s3(t)=fftshift(s3(t));
filter1=zeros(1,(4170/2)-417);
filter3=ones(1,2*417);
filter=[filter1 filter3 filter1] ;
g3_(t)=s3(t).*filter(t);

ecg3(t)=real(ifft(ifftshift(g3_(t)),'symmetric'));
%inverse fourier transform
figure(5)

plot(T(t),ecg3(t),'LineWidth',0.1);
title('signal after increasing SNR');

%% Asmaa's part
%-------- 4. Finding the heart rate using autocorrelation -----------------
AdultRate=90;
%ecg3
figure (6)
ecg3norm = ecg3 - mean(ecg3);
[autocor,lags] = xcorr(ecg3norm,ecg3norm);
autocor=autocor(round(length(autocor)/2):end);
lags=lags(round(length(lags)/2):end);

[peakShort,locationShort] = findpeaks(autocor);
short = mean(diff(locationShort))/fs;
[RRpeak,locationR] = findpeaks(autocor,'MinPeakDistance',ceil(short*fs),'MinPeakheight',6);
RRinterval = mean(diff(locationR))/(0.6*fs);

HeartRateECG3=AdultRate/RRinterval

plot(lags/fs,autocor);
hold on
plot(lags(locationShort)/fs,peakShort,'or', lags(locationR)/fs,RRpeak,'vk');
title('ECG3 autocorrelation')

%ecg2
figure (7)
ecg2norm = ecg2 - mean(ecg2);
[autocor,lags] = xcorr(ecg2norm,ecg2norm);
autocor=autocor(round(length(autocor)/2):end);
lags=lags(round(length(lags)/2):end);

[peakShort,locationShort] = findpeaks(autocor);
short = mean(diff(locationShort))/fs;
[RRpeak,locationR] = findpeaks(autocor,'MinPeakDistance',ceil(short*fs),'MinPeakheight',6);
RRinterval = mean(diff(locationR))/(0.6*fs);

HeartRateECG2=AdultRate/RRinterval

plot(lags/fs,autocor);
hold on
plot(lags(locationShort)/fs,peakShort,'or', lags(locationR)/fs,RRpeak,'vk');
title('ECG2 autocorrelation')

%ecg
figure (8)
ecgnorm = ecg - mean(ecg);
[autocor,lags] = xcorr(ecgnorm,ecgnorm);
autocor=autocor(round(length(autocor)/2):end);
lags=lags(round(length(lags)/2):end);

[peakShort,locationShort] = findpeaks(autocor);
short = mean(diff(locationShort))/fs;
[RRpeak,locationR] = findpeaks(autocor,'MinPeakDistance',ceil(short*fs),'MinPeakheight',6);
RRinterval = mean(diff(locationR))/(0.6*fs);
HeartRateECG=AdultRate/RRinterval


plot(lags/fs,autocor);
hold on
plot(lags(locationShort)/fs,peakShort,'or', lags(locationR)/fs,RRpeak,'vk');
title('ECG autocorrelation')

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% another solution %%%%%%%%%%%%%%%%%%%%%%%
load('ecg.mat');
t=1:1:length(ecg);
T=0:1/fs:length(ecg);
plot(T(t),ecg(t));

xlabel('milliseconds')
ylabel('millivolts')

% Remove trend from data
detrendedECG = detrend(ecg,5);
% Visualize results
clf
plot(ecg,'Color',[109 185 226]/255,'DisplayName','Input data')
hold on
plot(detrendedECG,'Color',[0 114 189]/255,'LineWidth',1.5,...
    'DisplayName','Detrended data')
plot(ecg-detrendedECG,'Color',[217 83 25]/255,'LineWidth',1,...
    'DisplayName','Trend')
hold off
legend
xlabel('milliseconds')
ylabel('millivolts')


%Find the R Wave Maximums
% Find local maxima
ismax = islocalmax(detrendedECG,'MinProminence',0.9);
% Visualize results
clf
plot(detrendedECG,'Color',[109 185 226]/255,'DisplayName','Input data')
hold on
% Plot local maxima
plot(find(ismax),detrendedECG(ismax),'^','Color',[217 83 25]/255,...
    'MarkerFaceColor',[217 83 25]/255,'DisplayName','Local maxima')
title(['Number of extrema: ' num2str(nnz(ismax))])
hold off
legend
xlabel('milliseconds')
ylabel('millivolts')

%Calculate the Heart Rate
maxIndices = find(ismax);
msPerBeat = mean(diff(maxIndices));
heartRate = 60*(1000/msPerBeat)
%}

%% conintuing of Asmaa's part
%-------- 5. Finding the QRS complex: -----------------

figure ; 
pan_tompkin(ecg, fs,1);

figure; 
pan_tompkin(ecg2, fs,1);

figure; 
pan_tompkin(ecg3, fs,1);

%{    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the detailed solution with steps %%%%%%%%%%%
if ~isvector(ecg)
  error('ecg must be a row or column vector');
end

gr = 1;   % on default the function always plots

ecg = ecg(:); % vectorize

% ======================= Initialize =============================== %
delay = 0;
skip = 0;                                                                  % becomes one when a T wave is detected
m_selected_RR = 0;
mean_RR = 0;
ser_back = 0; 
ax = zeros(1,6);

% ============ Noise cancelation(Filtering)( 5-15 Hz) =============== %%
if fs == 200
% ------------------ remove the mean of Signal -----------------------%
  ecg = ecg - mean(ecg);
% ==== Low Pass Filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2 ==== %%

   Wn = 12*2/fs;
   N = 3;                                                                  % order of 3 less processing
   [a,b] = butter(N,Wn,'low');                                             % bandpass filtering
   ecg_l = filtfilt(a,b,ecg); 
   ecg_l = ecg_l/ max(abs(ecg_l));
 % ======================= start figure ============================= %%
 
    figure;
    ax(1) = subplot(321);plot(ecg);axis tight;title('Raw signal');
    ax(2)=subplot(322);plot(ecg_l);axis tight;title('Low pass filtered');

% ==== High Pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1)) ==== %%
   Wn = 5*2/fs;
   N = 3;                                                                  % order of 3 less processing
   [a,b] = butter(N,Wn,'high');                                            % bandpass filtering
   ecg_h = filtfilt(a,b,ecg_l); 
   ecg_h = ecg_h/ max(abs(ecg_h));

    ax(3)=subplot(323);plot(ecg_h);axis tight;title('High Pass Filtered');

else
% bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
 f1=5;                                                                      % cuttoff low frequency to get rid of baseline wander
 f2=15;                                                                     % cuttoff frequency to discard high frequency noise
 Wn=[f1 f2]*2/fs;                                                           % cutt off based on fs
 N = 3;                                                                     % order of 3 less processing
 [a,b] = butter(N,Wn);                                                      % bandpass filtering
 ecg_h = filtfilt(a,b,ecg);
 ecg_h = ecg_h/ max( abs(ecg_h));
 
  ax(1) = subplot(3,2,[1 2]);plot(ecg);axis tight;title('Raw Signal');
  ax(3)=subplot(323);plot(ecg_h);axis tight;title('Band Pass Filtered');
 
end
% ==================== derivative filter ========================== %%
% ------ H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2)) --------- %
if fs ~= 200
 int_c = (5-1)/(fs*1/40);
 b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
 b = [1 2 0 -2 -1].*(1/8)*fs;   
end

 ecg_d = filtfilt(b,1,ecg_h);
 ecg_d = ecg_d/max(ecg_d);

  ax(4)=subplot(324);plot(ecg_d);
  axis tight;
  title('Filtered with the derivative filter');

% ========== Squaring nonlinearly enhance the dominant peaks ========== %%
 ecg_s = ecg_d.^2;

  ax(5)=subplot(325);
  plot(ecg_s);
  axis tight;
  title('Squared');


%============  Moving average ================== %%
%-------Y(nt) = (1/N)[x(nT-(N - 1)T)+ x(nT - (N - 2)T)+...+x(nT)]---------%
ecg_m = conv(ecg_s ,ones(1 ,round(0.150*fs))/round(0.150*fs));
delay = delay + round(0.150*fs)/2;


  ax(6)=subplot(326);plot(ecg_m);
  axis tight;
  title('Averaged with 30 samples length,Black noise,Green Adaptive Threshold,RED Sig Level,Red circles QRS adaptive threshold');
  axis tight;
%}



