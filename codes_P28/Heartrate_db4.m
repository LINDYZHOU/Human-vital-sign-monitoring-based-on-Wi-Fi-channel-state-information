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
frame_real_amp=framing(Rx1{1},sequence,frame_length,frame_interval);
num_f = length(frame_real_amp);
frame_real_ph=framing(Rx1{2},sequence,frame_length,frame_interval);
num_person=0;
multiperson_mode=2;
%% Each frame
for n= 1:num_f
    % Read each frame
    frame=frame_real_amp{n};
    %% people counting
    num_person_last=num_person;
    %last frame num_person value
    num_person=multi_person(frame_real_ph{n});
    % num_person=0: only one person
    % num_person=1: more than one person
    % num_person=2: not ax protocols
    if num_person < 2
        multiperson_mode = num_person_last & num_person;
        % continous 2 frames are determined there are multi-person under monitoring
        % the system should be carried on under multiperson mode
    end
    if multiperson_mode == 1
        error('是多人实验场景！')
    end
    %% SSNR-based subcarrier selection
    % x_b is for breath, x_h is for heartbeat.
    [csi_selected_h,selection_h{n}] = subcarrier_selection_ssnr_h(frame,10);
    [num_pkg,num_selected_h]=size(csi_selected_h);
    %% Bandpass filter for heart rate only
    csi_bpass_h=filter(bandpass_hb,csi_selected_h);
    % Time sequence before sampling (for heart rate).
    t_db4 = (0:0.1*num_pkg-1).*0.05;
    %% EMD, only for heartbeat
    % decompose heartbeat signal
%     emd_csi_h = emd_hb(csi_bpass_h,1);
    %% Wavelet Transform
    csi_wavelet_h_db4=[];
    for i = 1:num_selected_h
        csi_wavelet_h_db4(:,i) = wavelet_reference(csi_bpass_h(:,i));
    end
    %% Variance-based subcarrier group selection
    csi_var_h = VarianceAnalysis(csi_wavelet_h_db4,1);
%     figure();
%     plot(csi_var_h);
    %% Rate estimation
    % Heart rate. No PCA, each selected subcarrier is calculated.
    avgsmooth_h = [];
    for i = 1:size(csi_var_h,2)
        %Smooth
        avgsmooth_h(:,i)=smooth(csi_var_h(:,i),'rloess',50);
        avgsmooth_h(:,i)=avgsmooth_h(:,i)./(max(avgsmooth_h(:,i))-min(avgsmooth_h(:,i)));
        %Spectrum analysis based on DFT
        a_heart_rate_dft(n,i) = dft_estimation(avgsmooth_h(10:end,i),20);
        %Peaks and valleys interval analysis
        a_heart_rate_pp(n,i) = pp_estimation_db4(t_db4,avgsmooth_h(:,i));
    end
end
% Heart rate calculation.
% Remove the NaN values.
a_heart_rate_pp(a_heart_rate_pp==0)=NaN;
a_heart_rate_dft(a_heart_rate_dft==0)=NaN;
% The average value as result.
A_heart_rate_pp=mean(a_heart_rate_pp,"all",'omitnan');
A_heart_rate_dft = mean(a_heart_rate_dft,"all",'omitnan');
fprintf('单人！Heart rate: %fHz.\n', A_heart_rate_pp);