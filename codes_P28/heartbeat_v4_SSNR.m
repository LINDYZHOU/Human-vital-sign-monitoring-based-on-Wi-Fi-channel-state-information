% This project is done by Zhou Xiaolin.
% For heart rate monitoring, with SSNR.
% Email: 1477877598@qq.com

%% Read the CSI file
[Rx1,Rx2,sequence]=CSIdatapro(ax40_regular60s_zxl,'hesu');
% Extract the CSI data from RXs. Interpolate at the same time.

%% Framing
%every 10 seconds is a frame
frame_length=10;
%every 5 seconds produce a new frame
frame_interval=5;
%sampling rate is 200Hz
rate=200;

% Only the 1st Rx is used. Producing the frames.
frame_real=framing(Rx1{1},sequence,frame_length,frame_interval);
num_f = length(frame_real);
snr_last_selection=1:10;

for n= 1:num_f
    % Read each frame
    frame=frame_real{n};

    %     % Hampel filter
    %     frame = hampel(frame);

    %% SSNR-based & variance subcarrier selection
    %     [csi_selected_b,selection_b{n},location_b{n},snr_b{n}] = subcarrier_selection_ssnr_h(frame,10);
    variances = var(frame,0,1);
    [~,var_selected]=sortrows(variances','descend');
    csi_var_selected_h=frame(:,var_selected(1));

    [ssnr_csi_selected_h,selection_h{n}] = subcarrier_selection_ssnr_h(frame,10);
    [num_pkg,num_selected_h]=size(ssnr_csi_selected_h);

    %% SNR-based subcarrier selection
    [snr_csi_selected,snr_selection,~,~]=subcarrier_selection_pro(frame,0,10);
    if isempty(snr_selection)
        snr_csi_selected=frame(:,snr_last_selection);
    else
        snr_last_selection=snr_selection;
    end
    %     csi_filter=lowpass(snr_csi_selected,5,10,200);
    [snr_num_pkg,snr_num_selected]=size(snr_csi_selected);


    %% Bandpass filter for heart rate only
    csi_bpass_h=filter(bandpass_hb,ssnr_csi_selected_h);
    % x_h_bw means the heart signals from respiratory subcarrier selection.
    csi_bpass_h_var=filter(bandpass_hb,csi_var_selected_h);

    %     num_selected_h = size(csi_selected_h,2);

    % Time sequence before sampling (for heart rate).
    t = (0:num_pkg-1).*0.005;

    %% band-pass filter for breathing
    csiWavelet=[];
    for i = 1:snr_num_selected
        csiWavelet(:,i) = filter (filter5,snr_csi_selected(:,i));
    end
    [csi_sampled,t_sampling] = down_sample(csiWavelet,200,10);
    csiVar = VarianceAnalysis(csi_sampled,1);
    [coeff,score,latent] = pca(csiVar);
    Xcentered = score*coeff';
    avgsmooth=smooth(score(:,1),'rloess',50);

    respiratory_rate(n) = respirationrate_estimation(avgsmooth,t_sampling,frame_length);
    %% Downsampling, only for breath
    %     [csi_sampling,t_sampling] = down_sample(csi_selected_b,200,10);

    %% EMD, only for heartbeat
    % decompose heartbeat signal
    emd_csi_h = emd_hb(csi_bpass_h,1);
    emd_csi_h_var = emd_hb(csi_bpass_h_var,1);

    %% Wavelet Transform
    %     csi_wavelet_b=[];
    %     for i = 1:num_selected_b
    %         csi_wavelet_b(:,i) = wavelet_breathe (csi_sampling(:,i),6,'db6',location_b{n}(i));
    %     end

    csi_wavelet_h=[];
    for i = 1:num_selected_h
        csi_wavelet_h(:,i) = wavelet_heart(emd_csi_h(:,i),7,'coif4');
    end

        csi_wavelet_h_var=[];
            csi_wavelet_h_var = wavelet_heart(emd_csi_h_var,7,'coif4');

    %% Variance-based subcarrier group selection
    %     csi_var_b = VarianceAnalysis(csi_wavelet_b,1);
    csi_var_h = VarianceAnalysis(csi_wavelet_h,1);
    %     csi_var_h_bw = VarianceAnalysis(csi_wavelet_h_bw,1);

    %% Rate estimation
    % Heart rate. No PCA, each selected subcarrier is calculated.
    avgsmooth_h = [];
    for i = 1:size(csi_var_h,2)
        %Smooth
        avgsmooth_h(:,i)=smooth(csi_var_h(:,i),'rloess',50);
        avgsmooth_h(:,i)=avgsmooth_h(:,i)./(max(avgsmooth_h(:,i))-min(avgsmooth_h(:,i)));
        %Spectrum analysis based on DFT
        a_heart_rate_dft(n,i) = dft_estimation(avgsmooth_h(100:end,i),200);
        %Peaks and valleys interval analysis
        a_heart_rate_pp(n,i) = pp_estimation(t(100:end),avgsmooth_h(100:end,i),1);
    end

        avgsmooth_h_var = [];
            %Smooth
            avgsmooth_h_var=smooth(csi_wavelet_h_var,'rloess',50);
            avgsmooth_h_var=avgsmooth_h_var./(max(avgsmooth_h_var)-min(avgsmooth_h_var));
            %Spectrum analysis based on DFT
%             a_heart_rate_dft_bw(n,i) = dft_estimation(avgsmooth_h_var(100:end,i),200);
            %Peaks and valleys interval analysis
            a_heart_rate_pp_var(n) = pp_estimation(t(100:end),avgsmooth_h_var(100:end),1);

    %     % Respiratory rate.
    %     % PCA Principal Component Analysis (PCA).
    %     [~,score_b,~] = pca(csi_var_b);
    %     avgsmooth_b=smooth(score_b(:,1),'rloess',50);
    %     avgsmooth_b=avgsmooth_b./(max(avgsmooth_b)-min(avgsmooth_b));
    %     % Spectrum analysis based on DFT
    %     a_breath_rate_dft(n) = dft_estimation(avgsmooth_b,20);
    %     % Peaks and valleys interval analysis
    %     a_breath_rate_pp(n) = pp_estimation(t_sampling,avgsmooth_b,0);

end

% % Remove the outliers of respiratory rates from frames.
% a_breath_rate_dft_fd = rmoutliers(a_breath_rate_dft,'quartiles');
% a_breath_rate_pp_fd = rmoutliers(a_breath_rate_pp,'quartiles');
% % The average value as result.
% A_breath_rate_pp=mean(a_breath_rate_pp_fd,"all",'omitnan');
% A_breath_rate_dft = mean(a_breath_rate_dft_fd,"all",'omitnan');
% fprintf('Respiratory rate: DFT:  %fHz.   PEAK GAP: %fHz.\n', A_breath_rate_dft, A_breath_rate_pp);

% Heart rate calculation.
% Remove the NaN values.
a_heart_rate_pp(a_heart_rate_pp==0)=NaN;
a_heart_rate_dft(a_heart_rate_dft==0)=NaN;
% The average value as result.
A_heart_rate_pp=mean(a_heart_rate_pp,"all",'omitnan');
A_heart_rate_dft = mean(a_heart_rate_dft,"all",'omitnan');
fprintf('Heart rate: DFT: %fHz.PEAK: %fHz.\n', A_heart_rate_dft, A_heart_rate_pp);
A_respiratory=mean(respiratory_rate);
% fprintf('breathing rate via lowpass: %fHz',A_respiratory);

% % Remove the NaN values.
a_heart_rate_pp_var(a_heart_rate_pp_var==0)=NaN;
% a_heart_rate_dft_bw(a_heart_rate_dft_bw==0)=NaN;
% % The average value as result.
A_heart_rate_pp_var=mean(a_heart_rate_pp_var,"all",'omitnan');
% A_heart_rate_dft_bw = mean(a_heart_rate_dft_bw,"all",'omitnan');
fprintf('Heart rate largest var: PEAK GAP: %fHz.\n',  A_heart_rate_pp_var);


