%% read the CSI file
[Rx1,Rx2,sequence]=CSIdatapro(ax40_regular60s_wrl,'hesu');
% csi1 is the mag value of the file.
csi1 = [Rx1{1},Rx2{1}];
% [num_s,num_c]=size(csi1);%number of sample, number of subcarriers

%% framing
%every 10 seconds is a frame
frame_length=10;
%every 5 seconds produce a new frame
frame_interval=5;
%sampling rate is 200Hz
rate=200;

frame_real=framing(Rx1{1},sequence,frame_length,frame_interval);
% frame_real=framing(csi1,sequence,frame_length,frame_interval);
num_f = length(frame_real);
% num_f=round((num_s-frame_length*rate)/(rate*frame_interval));%number of frame

for n= 1:num_f
    %produce frame
    frame=frame_real{n};
    %% subcarrier selection
    [csi_selected_b,selection_b{n},location_b{n},snr_b{n}] = subcarrier_selection_pro(frame,0,10);
    % select better sensetive subcarriers from Spectrum
    [num_pkg,num_selected_b]=size(csi_selected_b);

    [csi_selected_h,selection_h{n},location_h{n},snr_h{n}] = subcarrier_selection_pro(frame,1,10);
    %selection for heartbeat

    %% bandpass filter
    csi_bpass_h=filter(bandpass_hb,csi_selected_b);
    csi_bpass_h2=filter(bandpass_hb,csi_selected_h);

    num_selected_h = size(csi_selected_h,2);
    
    t = (0:num_pkg-1).*0.005;
    % figure()
    % plot(t,csi_selected);
    % xlabel('Time')
    % ylabel('Amplitude')
    % title("After sub-selected")
    % grid on
    
    %% downsample, only for breath
    [csi_sampling,t_sampling] = down_sample(csi_selected_b,200,10);
%     figure();
%     plot(t,emd_csi_h);
%     xlabel('Time')
%     ylabel('Amplitude')
%     title('after EMD')
%     grid on
    
    %% EMD, only for heartbeat
    % decompose heartbeat signal
    emd_csi_h = emd_hb(csi_bpass_h,1);
    emd_csi_h2 = emd_hb(csi_bpass_h2,1);

    %% wavelet
    csi_wavelet_b=[];
    for i = 1:num_selected_b
        csi_wavelet_b(:,i) = wavelet_breathe (csi_sampling(:,i),6,'db6',location_b{n}(i));
    end
    
    csi_wavelet_h=[];
    for i = 1:num_selected_b
        csi_wavelet_h(:,i) = wavelet_heart(emd_csi_h(:,i),7,'coif4');
    end
    csi_wavelet_h2=[];
    for i = 1:num_selected_h
        csi_wavelet_h2(:,i) = wavelet_heart(emd_csi_h2(:,i),7,'coif4');
    end

    %% Variance
    csi_var_b = VarianceAnalysis(csi_wavelet_b,1);
    csi_var_h = VarianceAnalysis(csi_wavelet_h,1);
    csi_var_h2 = VarianceAnalysis(csi_wavelet_h2,1);
    
%     figure();
%     subplot(211)
%     plot(t_sampling,csi_var_b);
%     grid on
%     xlabel('Time')
%     ylabel('Amplitude')
% %     xlim([1,2]);
%     title("Before PCA")
%     subplot(212);
%     plot(t_sampling,avgsmooth_b); 
%     xlabel('Time');
%     ylabel('Amplitude');
% %     xlim([1,2]);
%     title('After PCA');
%     grid on
    
    %% PCA Principal Component Analysis (PCA) on raw data.
    [~,score_b,~] = pca(csi_var_b);
%     [~,score_h,~] = pca(csi_var_h);
%     [~,score_h2,~] = pca(csi_var_h2);
    % Xcentered = score*coeff';
%     figure();
%     plot(t_sampling,score_b(:,1));
%     xlabel('Time(s)')
%     ylabel('Amplitude')
%     title("After PCA")
%     grid on

    %% estimation
    %breath
%     avgsmooth_b = [];
%     for i = 1:size(csi_var_b,2)
%         %Smooth
%         avgsmooth_b(:,i)=smooth(csi_var_b(:,i),'rloess',50);
%         avgsmooth_b(:,i)=avgsmooth_b(:,i)./(max(avgsmooth_b(:,i))-min(avgsmooth_b(:,i)));
%         %DFT
%         a_breath_rate_dft(n,i) = dft_estimation(avgsmooth_b(10:end,i),20);
%         %PEAK-PEAK
%         a_breath_rate_pp(n,i) = pp_estimation(t_sampling(10:end),avgsmooth_b(10:end,i),0);
%     end
    %heart
    avgsmooth_h = [];
    for i = 1:size(csi_var_h,2)
        %Smooth
        avgsmooth_h(:,i)=smooth(csi_var_h(:,i),'rloess',50);
        avgsmooth_h(:,i)=avgsmooth_h(:,i)./(max(avgsmooth_h(:,i))-min(avgsmooth_h(:,i)));
        %DFT
        a_heart_rate_dft(n,i) = dft_estimation(avgsmooth_h(100:end,i),200);
        %PEAK-PEAK
        a_heart_rate_pp(n,i) = pp_estimation(t(100:end),avgsmooth_h(100:end,i),1);
    end
    %heart2
    avgsmooth_h2 =[];
    for i = 1:size(csi_var_h2,2)
        %Smooth
        avgsmooth_h2(:,i)=smooth(csi_var_h2(:,i),'rloess',50);
        avgsmooth_h2(:,i)=avgsmooth_h2(:,i)./(max(avgsmooth_h2(:,i))-min(avgsmooth_h2(:,i)));
        %DFT
        a_heart_rate_dft2(n,i) = dft_estimation(avgsmooth_h2(100:end,i),200);
        %PEAK-PEAK
        a_heart_rate_pp2(n,i) = pp_estimation(t(100:end),avgsmooth_h2(100:end,i),1);
    end

    avgsmooth_b=smooth(score_b(:,1),'rloess',50);
    avgsmooth_b=avgsmooth_b./(max(avgsmooth_b)-min(avgsmooth_b));
%     avgsmooth_h=smooth(score_h(:,1),'rloess',50);
%     avgsmooth_h=avgsmooth_h./(max(avgsmooth_h)-min(avgsmooth_h));
%     avgsmooth_h2=smooth(score_h2(:,1),'rloess',50);
%     avgsmooth_h2=avgsmooth_h2./(max(avgsmooth_h2)-min(avgsmooth_h2));
    
    %% DFT
    a_breath_rate_dft(n) = dft_estimation(avgsmooth_b,20);

%     a_heart_rate_dft(n) = dft_estimation(avgsmooth_h,200);
%     a_heart_rate_dft2(n) = dft_estimation(avgsmooth_h2,200);
    
    %% PeakPeak estimation
    a_breath_rate_pp(n) = pp_estimation(t_sampling,avgsmooth_b,0);

%     a_heart_rate_pp(n) = pp_estimation(t,avgsmooth_h,1);
%     a_heart_rate_pp2(n) = pp_estimation(t,avgsmooth_h2,1);
    
end

a_breath_rate_dft_fd = rmoutliers(a_breath_rate_dft);
a_breath_rate_pp_fd = rmoutliers(a_breath_rate_pp);
% 
% a_heart_rate_dft_fd = rmoutliers(a_heart_rate_dft);
% a_heart_rate_pp_fd = rmoutliers(a_heart_rate_pp);
% 
% a_heart_rate_dft_fd2 = rmoutliers(a_heart_rate_dft2);
% a_heart_rate_pp_fd2 = rmoutliers(a_heart_rate_pp2);

% conbine these 2 types of heart rate
% a_heart_rate_dft_cb = [a_heart_rate_dft,a_heart_rate_dft2];
% a_heart_rate_dft_cb_fd = rmoutliers(a_heart_rate_dft_cb,'mean');

% a_breath_rate_pp(a_breath_rate_pp==0)=NaN;
% a_breath_rate_dft(a_breath_rate_dft==0)=NaN;
% a_breath_rate_dft_f = mean(a_breath_rate_pp,2,'omitnan');
A_breath_rate_pp=mean(a_breath_rate_pp_fd,"all",'omitnan');
A_breath_rate_dft = mean(a_breath_rate_dft_fd,"all",'omitnan');
fprintf('Respiratory rate: DFT:  %fHz.   PEAK GAP: %fHz.\n', A_breath_rate_dft, A_breath_rate_pp);

a_heart_rate_pp(a_heart_rate_pp==0)=NaN;
a_heart_rate_dft(a_heart_rate_dft==0)=NaN;
A_heart_rate_pp=mean(a_heart_rate_pp,"all",'omitnan');
A_heart_rate_dft = mean(a_heart_rate_dft,"all",'omitnan');
fprintf('Heart rate: DFT: %fHz.     PEAK GAP: %fHz.\n', A_heart_rate_dft, A_heart_rate_pp);

a_heart_rate_pp2(a_heart_rate_pp2==0)=NaN;
a_heart_rate_dft2(a_heart_rate_dft2==0)=NaN;
A_heart_rate_pp2=mean(a_heart_rate_pp2,"all",'omitnan');
A_heart_rate_dft2 = mean(a_heart_rate_dft2,"all",'omitnan');
fprintf('Heart rate2: DFT: %fHz.    PEAK GAP: %fHz.\n', A_heart_rate_dft2, A_heart_rate_pp2);
