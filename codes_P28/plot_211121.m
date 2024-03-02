% read the CSI file
csi1 = CSIdata (ac20Zone2Exp1Breath7,'vht','Mag');
% csi1 is the mag value of the file.
%[num_t,num_subcarries] = size(csi1);

% use lowpass filter to remove the high-frequency noise.
csi_filter=lowpass(csi1,Wp,Ws,fk);
% Wp: pass cut-off frequence 
% Ws: block cut-off frequence 
% fk: sample rate(200)

% select better sensetive subcarriers from center Spectrum.
[csi_selected,selection]=subcarrier_selection(csi_filter,10);
[num_pkg,num_selected]=size(csi_selected);

% % change the col and row
% csi11 = csi1';

% plot amplitude
t = (0:num_pkg-1).*0.005;
figure()
for i = 1:num_selected
    plot(t,csi_selected(:,i));
    hold on    
end
hold off
xlabel('Time')
ylabel('Amplitude')
title("After sub-selected")

% downsample
[csi_sampling,t_sampling] = down_sample(csi_selected,ori_srate,sample_times);
figure();
plot(t_sampling,csi_sampling);
xlabel('Time')
ylabel('Amplitude')
title('after Downsample')
axis([0,20,-20,20]);

% % wavelet, use'waveletAnalyzer'
% % VAR_1 is db6
% % VAR_2 is sym6
% subcarrier_4 = my_VAR_1(5,:)+my_VAR_1(6,:);
% figure()
% plot(t_sampling,csiWavelet);
% xlabel('Time')
% ylabel('Amplitude')

%wavelet
for i = 1:num_selected
    csiWavelet(:,i) = wavelet_breathe (csi_sampling(:,i),6,'db6');
end
figure()
plot(t_sampling,csiWavelet);
xlabel('Time')
ylabel('Amplitude')
title("After Wavelet")

% Variance
csiVar = VarianceAnalysis(csiWavelet,1);
figure();
plot(t_sampling,csiVar);
xlabel('Time')
ylabel('Amplitude')
title("After Variance selection")

% PCA Principal Component Analysis (PCA) on raw data.
[coeff,score,latent] = pca(csiVar);
Xcentered = score*coeff'
figure();
plot(t_sampling,Xcentered(:,1));
xlabel('Time')
ylabel('Amplitude')
title("After PCA")

% Smooth
avgsmooth=smooth(Xcentered(:,1),'rloess',50);

% respiratory rate
[pks, locs] = findpeaks(avgsmooth,t_sampling,'MinPeakWidth',0.5);
breathtimes=length(pks);
peakInterval = diff(locs);
respiratory_duration = mean(peakInterval);
respiratory_rate = respiratory_duration/20;
str = ["respiratory interval: (seconds)",num2str(peakInterval),"respiratory duration = (seconds)", num2str(respiratory_duration),"respiratory rate = (Hz)", num2str(respiratory_rate)];

figure();
plot(t_sampling,avgsmooth);
grid on
xlabel('Time')
ylabel('Amplitude')
title("After Smooth")
text(2,2,str);


