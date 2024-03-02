% read the CSI file
csi1 = CSIdata (ac20_hold10s_zxl,'vht','Mag');
% csi1 is the mag value of the file.
%[num_t,num_subcarries] = size(csi1);

csi_filter=lowpass(csi1,5,10,200);
% use lowpass filter to remove the high-frequency noise.
% Wp: pass cut-off frequence 
% Ws: block cut-off frequence 
% fk: sample rate(200)
% csi_filter = lowpass(input,Ws1,Wp1,Wp2,Ws2,fk)

[csi_selected,selection]=subcarrier_selection(csi_filter,10);
% select better sensetive subcarriers from Spectrum
[num_pkg,num_selected]=size(csi_selected);

% plot amplitude
t = (0:num_pkg-1).*0.005;
figure()
plot(t,csi_selected);
xlabel('Time')
ylabel('Amplitude')
title("After sub-selected")
grid on

% downsample
[csi_sampling,t_sampling] = down_sample(csi_selected,200,10);
figure();
plot(t_sampling,csi_sampling);
xlabel('Time')
ylabel('Amplitude')
title('after Downsample')
grid on

%wavelet
for i = 1:num_selected
    csiWavelet(:,i) = wavelet_breathe (csi_sampling(:,i),6,'db6');
end
figure()
plot(t_sampling,csiWavelet);
xlabel('Time')
ylabel('Amplitude')
title("After Wavelet")
grid on

% Variance
csiVar = VarianceAnalysis(csiWavelet,1);
figure();
plot(t_sampling,csiVar);
xlabel('Time')
ylabel('Amplitude')
title("After Variance selection")
grid on

% PCA Principal Component Analysis (PCA) on raw data.
[coeff,score,latent] = pca(csiVar);
% Xcentered = score*coeff';
figure();
plot(t_sampling,score(:,1));
xlabel('Time')
ylabel('Amplitude')
title("After PCA")
grid on

% Smooth
avgsmooth=smooth(score(:,1),'rloess',50);

% respiratory rate
avgsmooth=avgsmooth./(max(avgsmooth)-min(avgsmooth));
[~,~,~,hight] = findpeaks(avgsmooth,'MinPeakWidth',10);
[pks, locs] = findpeaks(avgsmooth,t_sampling,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
breathtimes=length(pks);
peakInterval = diff(locs);
respiratory_duration = mean(peakInterval);
respiratory_rate = respiratory_duration/20;
str = {"breathtimes = ",num2str(breathtimes),"respiratory interval: (seconds)",num2str(peakInterval),"respiratory duration = (seconds)", num2str(respiratory_duration),"respiratory rate = (Hz)", num2str(respiratory_rate)};

figure();
% findpeaks(avgsmooth,t_sampling,'MinPeakProminence',2,'Annotate','extents');
plot(t_sampling,avgsmooth);
text(1,1,str)
grid on
xlabel('Time')
ylabel('Amplitude')
title("After Smooth")



