%%  read the CSI file
csi1 = CSIdata_in (ac20_regular60s_hyh,'vht','Mag');
% csi1 is the mag value of the file.
%[num_t,num_subcarries] = size(csi1);
%sampling rate is 200Hz
rate=200;

[num_s,num_subcarriers]=size(csi1);%number of sample, number of subcarriers

%every 10 seconds is a frame
frame_length=10;
%every 5 seconds produce a new frame
frame_interval=5;
num_f=round((num_s-frame_length*rate)/(rate*frame_interval));%number of frame

%% produce frame
for n= 1:(num_f+1)
    if n ~= (num_f+1)
        frame=csi1((n-1)*(frame_interval*rate)+1:(n-1)*(frame_interval*rate)+frame_length*rate,:);
    else
        frame=csi1((n-1)*(frame_interval*rate)+1:end,:);
    end

[selected_h,selection{n},location{n},snr{n}] = subcarrier_selection_pro(frame,1,10);
% [selected_h,selection] = subcarrier_selection_h(rate,csi_filter,10);
%selection for heartbeat

csi_filter=filter(bandpass_hb,selected_h);

%     figure();     
% for csi_i=1:10
% N = 2^(2+nextpow2(length(csi_bpass_h)));
%     fft_sub=fft(csi_bpass_h(:,csi_i),N);
%     fftshift_sub=fftshift(fft_sub); 
%     f=(-N/2:N/2-1).*200./N;
%     plot(f,(abs(fftshift_sub).^2));
%     hold on;
% end
%     xlim([0,4]);
%     xlabel('频率/ Hz');
%     ylabel('功率');
%     title('带通滤波后的已选子载波');
%     grid on;  


[num_pkg,num_selected]=size(selected_h);

% plot amplitude
ts = (0:num_pkg-1).*0.005;
% figure()
% plot(t,csi_selected);
% xlabel('Time')
% ylabel('Amplitude')
% title("After sub-selected")
% grid on

% downsample
%[csi_sampling,t_sampling] = down_sample(csi_selected,200,10);
% figure();
% plot(t_sampling,csi_sampling);
% xlabel('Time')
% ylabel('Amplitude')
% title('after Downsample')
% grid on

%% EMD
% decompose heartbeat signal
emd_csi_h = emd_hb(csi_filter,1);

% figure();
%     subplot(211);
%         plot(csi_filter);
%         xlabel('time/20ms');
%         ylabel('amplitude');
%         title('filtered heartbeat signal');
%     subplot(212);
%         plot(emd_csi_h);
%         xlabel('time');
%         ylabel('amplitude');
%         title('filtered heartbeat signal by EMD');

%%  wavelet
wavelet_csi_h=[];
for i = 1:num_selected
    wavelet_csi_h(:,i) = wavelet_heart(emd_csi_h(:,i),7,'coif4');
end

figure();  
    subplot(211);
        plot(t,csi_selected_h);    
        grid on
        xlabel('Time');
        ylabel('Amplitude');
        title('Selected signal');
    subplot(212);
        plot(t,csi_wavelet_h); 
        grid on
        xlabel('time/100ms');
        ylabel('Amplitude');
        title('filtered heartbeat signal by WT');

%lowpass filter
% csiWavelet=[];
% for i = 1:num_selected
%      csiWavelet(:,i) = filter (filter5,csi_selected(:,i));
% end
% [csiWavelet,t_sampling] = down_sample(csiWavelet,200,10);
% figure()
% plot(t_sampling,csiWavelet);
% xlabel('Time')
% ylabel('Amplitude')
% title("After Wavelet")
% grid on

% Variance
csiVar = VarianceAnalysis(wavelet_csi_h,1);
% figure();
% plot(ts,csiVar);
% xlabel('Time')
% ylabel('Amplitude')
% title("After Variance selection")
% grid on

% PCA Principal Component Analysis (PCA) on raw data.
% [~,score,~] = pca(csiVar);
% Xcentered = score*coeff';
% figure();
% plot(ts,score(:,1));
% xlabel('Time(s)')
% ylabel('Amplitude')
% title("After PCA")
% grid on

for i = 1:size(csiVar,2)
% Smooth
avgsmooth(:,i)=smooth(csiVar(:,i),'rloess',50);
avgsmooth(:,i)=avgsmooth(:,i)./(max(avgsmooth(:,i))-min(avgsmooth(:,i)));
%% DFT
a_heart_rate_dft(n,i) = dft_estimation(avgsmooth(100:end,i),200);
a_heartrate_pp(n,i) = pp_estimation(ts(100:end),avgsmooth(100:end,i),1);
end


% N = 2^(4+nextpow2(num_pkg));
% fft_sub=fft(avgsmooth,N);
% fftshift_sub=fftshift(fft_sub); 
% f=(-N/2:N/2-1).*rate./N;
% psd = (abs(fftshift_sub(1:N)).^2)./N;
% figure();     
% plot(f,(abs(fftshift_sub).^2));
% hold on
% plot(f,psd);
% xlim([-5,5]);
% xlabel('frequncy / Hz');
% ylabel('amplitude');
% title('DFT of heartbeat signal pca');
% grid on;  
% [~,index_h]=max(psd);
% a_Heart_rate_dft(n) =abs(f(index_h));   

% heartbeat rate estimation
% [~,~,~,hight] = findpeaks(avgsmooth,'MinPeakWidth',10);
% [pks_peak, locs_peak] = findpeaks(avgsmooth,ts,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
% [pks_valley, locs_valley] = findpeaks(-1.*avgsmooth,ts,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
% %if there is 2 breaths in 1 frame
% if length(pks_peak) >= 2 && length(pks_valley) >= 2
%     % no less than 2 times
%     %combine location of peaks and valleys
%     locs=[locs_peak(1:min([length(locs_peak),length(locs_valley)]));locs_valley(1:min([length(locs_peak),length(locs_valley)]))];
%     %find the difference time of the adjacent peak and valley, which is
%     %half of the duration
%     Interval = abs(diff(locs,1,1));
%     respiratory_duration = mean(Interval);
%     respiratory_rate(n) = 0.5/respiratory_duration;
% else
%     % less than 2 times
%     % mirror the avgsmmoth twice  
%     breath_wave=[fliplr(avgsmooth')';avgsmooth;fliplr(avgsmooth')'];
% %     figure();
% %     plot(breath_wave);
%     % 3 times longer and then findpeaks
%     % make man-made duration to hlep findpeaks
%     [~,~,~,hight] = findpeaks(breath_wave,'MinPeakWidth',10);
%     pks = findpeaks(breath_wave,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
%     vals = findpeaks(-1.*breath_wave,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
%     breath_times=(length(pks)+length(vals))/6;
%     respiratory_rate(n)=breath_times/frame_length;
% end

% text(6,0,['respiratory rate:',num2str(respiratory_rate(n))])
%str = {"breathtimes = ",num2str(breathtimes(i)),"respiratory interval: (seconds)",num2str(peakInterval),"respiratory duration = (seconds)", num2str(respiratory_duration),"respiratory rate = (Hz)", num2str(respiratory_rate(i))};

% figure();
% % findpeaks(avgsmooth,t_sampling,'MinPeakProminence',2,'Annotate','extents');
% plot(t_sampling,avgsmooth);
% % text(1,1,str)
% % grid on
% xlabel('Time')
% ylabel('Amplitude')
% title(["Frame ";num2str(n)])
end

A_heartrate_pp=mean(a_heartrate_pp,"all");
fprintf('PEAK GAP: Heart rate: %fHz. \n', A_heartrate_pp);
A_Heart_rate_dft = mean(a_heart_rate_dft,"all");
fprintf('DFT: Heart rate: %fHz. \n', A_Heart_rate_dft);

