% read the CSI file
csi1 = CSIdata_in (ac40_regular60s_wrl,'vht','Mag');
% csi1 is the mag value of the file.
[num_s,num_c]=size(csi1);%number of sample, number of subcarriers


%% framing
%every 10 seconds is a frame
frame_length=10;
%every 5 seconds produce a new frame
frame_interval=5;
%sampling rate is 200Hz
rate=200;
num_f=round((num_s-frame_length*rate)/(rate*frame_interval));%number of frame

for n= 1:(num_f+1)
    %produce frame
    if n ~= (num_f+1)
        frame=csi1((n-1)*(frame_interval*rate)+1:(n-1)*(frame_interval*rate)+frame_length*rate,:);
    else
        frame=csi1((n-1)*(frame_interval*rate)+1:end,:);
    end

    % lowpass filter
%     csi_filter=filter(lowpass,frame);

    %% subcarrier selection
%     [csi_selected,selection,location,th]=subcarrier_selection_b(csi_filter,200);
    [csi_selected,selection(n,:),location(n,:),snr(n,:)] = subcarrier_selection_pro(frame,10);
% for csi_i=1:4
% N = 2^(4+nextpow2(length(csi_selected)));
%     fft_sub=fft(csi_selected(:,csi_i),N);
%     fftshift_sub=fftshift(fft_sub); 
%     f=(-N/2:N/2-1).*200./N;
%         figure();     
%     plot(f,(abs(fftshift_sub).^2));
% %     % plot(f,psd);
%     xlim([-5,5]);
%     xlabel('frequncy / Hz');
%     ylabel('amplitude');
%     title('DFT of heartbeat signal pca');
%     grid on;  
% end


    % select better sensetive subcarriers from Spectrum
    [num_pkg,num_selected]=size(csi_selected);
%     fprintf('NO. %d frame. Selection: \n', n);
%     selection(n,:)

    t = (0:num_pkg-1).*0.005;

%     figure()
%     subplot(211)
%     plot(t,csi_wavelet_h2);
%     xlabel('时间/10s')
%     ylabel('幅度')
%     title("使用db4小波基")
%     grid on
%     subplot(212)
%     plot(t,csi_wavelet_h);
%     xlabel('时间/10s')
%     ylabel('幅度')
%     title("使用coif4小波基")
%     grid on
%     
    %% downsample, only for breath
    [csi_sampling,t_sampling] = down_sample(csi_selected,200,10);
    % figure();
    % plot(t_sampling,csi_sampling);
    % xlabel('Time')
    % ylabel('Amplitude')
    % title('after Downsample')
    % grid on
    
    %% wavelet
    csi_wavelet=[];
    for i = 1:num_selected
        csi_wavelet(:,i) = wavelet_breathe (csi_sampling(:,i),6,'db6',location(n,i));
        % approximation or detail coefficients
    end
    
    
    %% Variance
    csi_var = VarianceAnalysis(csi_wavelet,1);
    
%     figure();
%     plot(t_sampling,csi_var);
%     xlabel('Time')
%     ylabel('Amplitude')
%     title("After Variance selection")
%     grid on
    
    %% PCA Principal Component Analysis (PCA) on raw data.
    [~,score,~] = pca(csi_var);
    % Xcentered = score*coeff';
%     figure();
%     plot(t_sampling,avgsmooth_b);
%     xlabel('Time(s)')
%     ylabel('Amplitude')
%     title("After PCA")
%     grid on
    
    % Smooth
    avgsmooth=smooth(score(:,1),'rloess',50);
    avgsmooth=avgsmooth./(max(avgsmooth)-min(avgsmooth));
    
    %% DFT
%     N = 2^(4+nextpow2(length(t_sampling)));
%     fft_sub=fft(avgsmooth,N);
%     fftshift_sub=fftshift(fft_sub); 
%     f=(-N/2:N/2-1).*(rate/10)./N;
    %     figure();     
%     plot(f,(abs(fftshift_sub).^2));
%     % plot(f,psd);
%     xlim([-5,5]);
%     xlabel('frequncy / Hz');
%     ylabel('amplitude');
%     title('DFT of heartbeat signal pca');
%     grid on;  
%     psd = (abs(fftshift_sub(1:N)).^2)./N;
%     [~,index_h]=max(psd);
%     [~,index_h]=max((abs(fftshift_sub).^2));
%     a_breath_rate_dft(n) =abs(f(index_h));   
    a_breathrate_dft(n) = dft_estimation(avgsmooth,20);
    
    % heartbeat rate estimation
    a_breathrate_pp(n) = pp_estimation(t_sampling,avgsmooth,0);
    
end

a_breathrate_dft_fd = rmoutliers(a_breathrate_dft);
a_breathrate_pp_fd = rmoutliers(a_breathrate_pp);

A_breathrate_pp=mean(a_breathrate_pp_fd);
fprintf('PEAK GAP: Respiratory rate: %fHz. \n', A_breathrate_pp);
A_breathrate_dft = mean(a_breathrate_dft_fd);
fprintf('DFT: Respiratory rate: %fHz. \n', A_breathrate_dft);

% possibility of 2 person
% person = var(selection,0,2);
