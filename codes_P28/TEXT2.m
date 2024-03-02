% read the CSI file
csi1 = CSIdata (ax40Exp2MaleHoldBreath3,'hesu','Mag');
% csi1 is the mag value of the file.
%[num_t,num_subcarries] = size(csi1);

%every 10 seconds is a frame
frame_length=10;
%every 5 seconds produce a new frame
frame_interval=5;
%sampling rate is 200Hz
rate=200;
[num_s,num_c]=size(csi1);%number of sample, number of subcarriers
num_f=round((num_s-frame_length*rate)/(rate*frame_interval));%number of frame
for n= 1:(num_f+1)
    %produce frame
    if n ~= (num_f+1)
        frame=csi1((n-1)*(frame_interval*rate)+1:(n-1)*(frame_interval*rate)+frame_length*rate,:);
    else
        frame=csi1((n-1)*(frame_interval*rate)+1:end,:);
    end

csi_filter=lowpass(frame,5,10,200);
% use lowpass filter to remove the high-frequency noise.
% Wp: pass cut-off frequence 
% Ws: block cut-off frequence 
% fk: sample rate(200)
% csi_filter = lowpass(input,Ws1,Wp1,Wp2,Ws2,fk)

[csi_selected,selection,location]=subcarrier_selection(csi_filter,10);
% select better sensetive subcarriers from Spectrum
[num_pkg,num_selected]=size(csi_selected);

% plot amplitude
t = (0:num_pkg-1).*0.005;
% figure()
% plot(t,csi_selected);
% xlabel('Time')
% ylabel('Amplitude')
% title("After sub-selected")
% grid on

% downsample
[csi_sampling,t_sampling] = down_sample(csi_selected,200,10);
% figure();
% plot(t_sampling,csi_sampling);
% xlabel('Time')
% ylabel('Amplitude')
% title('after Downsample')
% grid on

% wavelet
csiWavelet=[];
for i = 1:num_selected
    csiWavelet(:,i) = wavelet_breathe (csi_sampling(:,i),6,'db6',location(i));
end

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
csiVar = VarianceAnalysis(csiWavelet,1);
% figure();
% plot(t_sampling,csiVar);
% xlabel('Time')
% ylabel('Amplitude')
% title("After Variance selection")
% grid on

% PCA Principal Component Analysis (PCA) on raw data.
[coeff,score,latent] = pca(csiVar);
Xcentered = score*coeff';
% figure();
% plot(t_sampling,score(:,1));
% xlabel('Time(s)')
% ylabel('Amplitude')
% % title("After PCA")
% grid on

% Smooth
avgsmooth=smooth(score(:,1),'rloess',50);

% respiratory rate
avgsmooth=avgsmooth./(max(avgsmooth)-min(avgsmooth));
[~,~,~,hight] = findpeaks(avgsmooth,'MinPeakWidth',10);
[pks_peak, locs_peak] = findpeaks(avgsmooth,t_sampling,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
[pks_valley, locs_valley] = findpeaks(-1.*avgsmooth,t_sampling,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
%if there is 2 breaths in 1 frame
if length(pks_peak) >= 2 && length(pks_valley) >= 2
    % no less than 2 times
    %combine location of peaks and valleys
    locs=[locs_peak(1:min([length(locs_peak),length(locs_valley)]));locs_valley(1:min([length(locs_peak),length(locs_valley)]))];
    %find the difference time of the adjacent peak and valley, which is
    %half of the duration
    Interval = abs(diff(locs,1,1));
    respiratory_duration = mean(Interval);
    respiratory_rate(n) = 0.5/respiratory_duration;
else
    % less than 2 times
    % mirror the avgsmmoth twice  
    breath_wave=[fliplr(avgsmooth')';avgsmooth;fliplr(avgsmooth')'];
%     figure();
%     plot(breath_wave);
    % 3 times longer and then findpeaks
    % make man-made duration to hlep findpeaks
    [~,~,~,hight] = findpeaks(breath_wave,'MinPeakWidth',10);
    pks = findpeaks(breath_wave,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
    vals = findpeaks(-1.*breath_wave,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
    breath_times=(length(pks)+length(vals))/6;
    respiratory_rate(n)=breath_times/frame_length;
end
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
A_respiratory=mean(respiratory_rate);


