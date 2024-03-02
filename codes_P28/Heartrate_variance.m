% This project is done by Zhou Xiaolin.
% For heart rate monitoring, with SSNR.
% Email: 1477877598@qq.com
%% Read the CSI file
[Rx1,~,sequence]=CSIdatapro(ax40_regular60s_zxl,'hesu');
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
% snr_last_selection=1:10;
%% Each frame
for n= 1:num_f
    % Read each frame
    frame=frame_real{n};
    %% variance subcarrier selection
    variances = var(frame,0,1);
    [~,var_selected]=sortrows(variances','descend');
    csi_var_selected_h=frame(:,var_selected(1));
    [num_pkg,num_selected_h]=size(csi_var_selected_h);
    %% Bandpass filter for heart rate only
    csi_bpass_h_var=filter(bandpass_hb,csi_var_selected_h);
    % Time sequence before sampling (for heart rate).
    t = (0:num_pkg-1).*0.005;
    %% EMD, only for heartbeat
    % decompose heartbeat signal
    emd_csi_h_var = emd_hb(csi_bpass_h_var,1);
    %% Wavelet Transform
    csi_wavelet_h_var=[];
    csi_wavelet_h_var = wavelet_heart(emd_csi_h_var,7,'coif4');
%     %% Variance-based subcarrier group selection
%     csi_var_h = VarianceAnalysis(csi_wavelet_h,1);
    %% Rate estimation
    % Heart rate. No PCA, each selected subcarrier is calculated.
    avgsmooth_h_var = [];
    %Smooth
    avgsmooth_h_var=smooth(csi_wavelet_h_var,'rloess',50);
    avgsmooth_h_var=avgsmooth_h_var./(max(avgsmooth_h_var)-min(avgsmooth_h_var));
    %Peaks and valleys interval analysis
    a_heart_rate_pp_var(n) = pp_estimation(t(100:end),avgsmooth_h_var(100:end),1);
end
%% Heart rate calculation.
% % Remove the NaN values.
a_heart_rate_pp_var(a_heart_rate_pp_var==0)=NaN;
% a_heart_rate_dft_bw(a_heart_rate_dft_bw==0)=NaN;
% % The average value as result.
A_heart_rate_pp_var=mean(a_heart_rate_pp_var,"all",'omitnan');
% A_heart_rate_dft_bw = mean(a_heart_rate_dft_bw,"all",'omitnan');
fprintf('Heart rate the largest var: %fHz.\n',  A_heart_rate_pp_var);