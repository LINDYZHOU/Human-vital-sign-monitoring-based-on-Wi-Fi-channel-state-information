% load csi data
[Rx1,Rx2,sequence]=CSIdatapro(ax40_regular60s_hyh,'hesu');
frame_length=10;
%every 5 seconds produce a new frame
frame_interval=5;
%sampling rate is 200Hz
rate=200;
frame1=framing(Rx1{1},sequence,frame_length,frame_interval);
l=length(frame1);
selection=[];

for n=1:l
    frame=frame1{n};
    selection_last=selection;
    [sub_selected,selection,location]=subcarrier_selection_pro(frame,10);
    if isempty(selection)
        selection=selection_last;
        sub_selected=frame(:,selection);
    end
    [num_pkg,num_selected]=size(sub_selected);
%     Fs=200;
%     for i=1:length(selection)
%         sub=sub_selected(:,i);
%         num_pkg=length(sub);
%         N=2^(nextpow2(num_pkg)+3);
%         f=Fs*(-N/2:(N/2-1))/N;
%         fft_sub=fft(sub,N);
%         fftshift_sub=fftshift(fft_sub);
%         plot(f,abs(fftshift_sub));
%         axis([-1,1,-inf,inf]);
%     end
[phase_sampling,t_sampling] = down_sample(sub_selected,200,10);

% csiWavelet=[];
% for i = 1:num_selected
%     csiWavelet(:,i) = wavelet_breathe (phase_sampling(:,i),6,'db6',location(i));
% end
csiWavelet_bre = wavelet_breathe_pro(phase_sampling,'db4',1);

% figure()
% plot(t_sampling,csiWavelet);
% xlabel('Time')
% ylabel('Amplitude')
% title("After Wavelet")
% grid on

csiVar_bre = VarianceAnalysis(csiWavelet_bre,1);

[coeff,score,latent] = pca(csiVar_bre);
Xcentered = score*coeff';

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
% figure();
% % findpeaks(avgsmooth,t_sampling,'MinPeakProminence',2,'Annotate','extents');
% plot(t_sampling,avgsmooth);

csiWavelet_heart = wavelet_breathe_pro(phase_sampling,'db4',2);
csiVar_heart = VarianceAnalysis(csiWavelet_heart,1);
[avgheart(n),heavist_heart(n),~]=heartrate_estimation(csiVar_heart);

end
A_respiratory=mean(respiratory_rate);
A_heartrate=[mean(avgheart),mean(heavist_heart)];
% %lowpass filter
% % Rx1_amp_filtered=filter(filter5,Rx1_amp);
% Rx1_amp_filtered=lowpass(Rx1_amp,5,10,200);