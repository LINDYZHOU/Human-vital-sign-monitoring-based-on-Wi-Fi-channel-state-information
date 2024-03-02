%% test the numerous subcarriers could help distinguish 2 persons situation
% or not
% load csi data
[Rx1,Rx2,sequence]=CSIdatapro(ax20_regular60s_65_12,'hesu');

% frame to reduce the inputting data scale
%every 10 seconds is a frame
frame_length=10;
%every 5 seconds produce a new frame
frame_interval=5;
%sampling rate is 200Hz
rate=200;
[num_s,num_c]=size(Rx1{1});%number of sample, number of subcarriers
num_f=round((num_s-frame_length*rate)/(rate*frame_interval));%number of frame
sta_dev=zeros((num_f+1),4);
for n= 1:(num_f+1)
    %produce frame
    if n ~= (num_f+1)
        frame_rx1_amp=Rx1{1}((n-1)*(frame_interval*rate)+1:(n-1)*(frame_interval*rate)+frame_length*rate,:);
        frame_rx1_pha=Rx1{2}((n-1)*(frame_interval*rate)+1:(n-1)*(frame_interval*rate)+frame_length*rate,:);
        frame_rx2_amp=Rx2{1}((n-1)*(frame_interval*rate)+1:(n-1)*(frame_interval*rate)+frame_length*rate,:);
        frame_rx2_pha=Rx2{2}((n-1)*(frame_interval*rate)+1:(n-1)*(frame_interval*rate)+frame_length*rate,:);
    else
        frame_rx1_amp=Rx1{1}((n-1)*(frame_interval*rate)+1:end,:);
        frame_rx1_pha=Rx1{2}((n-1)*(frame_interval*rate)+1:end,:);
        frame_rx2_amp=Rx2{1}((n-1)*(frame_interval*rate)+1:end,:);
        frame_rx2_pha=Rx2{2}((n-1)*(frame_interval*rate)+1:end,:);
    end

%calculate the 
% R_amp = corrcoef(frame_rx1_amp);
% R_pha = corrcoef(frame_rx1_pha);

num_sub=20;
%     use subcarrier selection
%     based on snr
%     figure();
    [csi_snr_amp,selection_amp_snr,~,snr_amp]=subcarrier_selection_pro(frame_rx1_amp,num_sub);
    R_snr_amp = corrcoef(csi_snr_amp);
%     subplot(2,3,1);
%     scatter(selection_amp_snr,snr_amp);
%     title('snr amp');
    [csi_snr_pha,selection_pha_snr,~,snr_pha]=subcarrier_selection_pro(frame_rx1_pha,num_sub);
    R_snr_pha = corrcoef(csi_snr_pha);
%     subplot(2,3,2);
%     scatter(selection_pha_snr,snr_pha);
%     title('snr pha');
    [csi_snr_phadiff,selection_phadiff_snr,~,snr_phadiff]=subcarrier_selection_pro(frame_rx1_pha-frame_rx2_pha,num_sub);
    R_snr_phadiff = corrcoef(csi_snr_phadiff);
%     subplot(2,3,3);
%     scatter(selection_phadiff_snr,snr_phadiff);
%     title('snr phadiff');
    % based on ssnr
    [csi_ssnr_amp,selection_amp_ssnr,~,ssnr_amp]=subcarrier_selection_ssnr(frame_rx1_amp,num_sub);
    R_ssnr_amp = corrcoef(csi_ssnr_amp);
%     subplot(2,3,4);
%     scatter(selection_amp_ssnr,ssnr_amp);
%     title('ssnr amp');
    [csi_ssnr_pha,selection_pha_ssnr,~,ssnr_pha]=subcarrier_selection_ssnr(frame_rx1_pha,num_sub);
    R_ssnr_pha = corrcoef(csi_ssnr_pha);
%     subplot(2,3,5);
%     scatter(selection_pha_ssnr,ssnr_pha);
%     title('ssnr pha');
    [csi_ssnr_phadiff,selection_phadiff_ssnr,~,ssnr_phadiff]=subcarrier_selection_ssnr(frame_rx1_pha-frame_rx2_pha,num_sub);
    R_ssnr_phadiff = corrcoef(csi_ssnr_phadiff);
%     subplot(2,3,6);
%     scatter(selection_phadiff_ssnr,ssnr_phadiff);
%     title('ssnr phadiff');

    % calculate stadard deviation
   sta_dev(n,:)=[std(selection_amp_snr),std(selection_pha_snr),std(selection_amp_ssnr),std(selection_pha_ssnr)];
%    R_min(n,:)=[min(R_snr_amp,[],"all"),min(R_snr_pha,[],"all"),min(R_snr_phadiff,[],"all"),min(R_ssnr_amp,[],"all"),min(R_ssnr_pha,[],"all"),min(R_ssnr_phadiff,[],"all")];
   R_mean(n,:)=[mean(R_snr_amp,'all'),mean(R_snr_pha,'all'),mean(R_snr_phadiff,'all'),mean(R_ssnr_amp,'all'),mean(R_ssnr_pha,'all'),mean(R_ssnr_phadiff,'all')];

end

R_ssnr=load(['R_mean','.mat']).R_ssnr;
R_ssnr=[R_ssnr;R_mean(:,5)];
save('R_mean','R_ssnr');

