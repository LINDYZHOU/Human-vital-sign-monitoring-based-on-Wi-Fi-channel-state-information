%% Read the CSI data from .csi file

csi_raw = CSIdata (ax40_regular60s_wrl,'hesu','Mag');
Fs_h = 200;
[num_pkg,num_subcarriers]=size(csi_raw);
ts = (1:num_pkg) .*0.005;

%% Hampel filter  
% Remove outliners and replace them with local median    

% filter heartbeat signal    
[hamed_csi_h, outlier_h]=hampel(csi_raw,3);     %the length of window of both side is 3

% csi_filter=lowpass(csi_raw,2,5,200);
csi_filter=filter(lowpass_hb,csi_raw);
% use lowpass filter to remove the high-frequency noise.
% Wp: pass cut-off frequence 
% Ws: block cut-off frequence 
% fk: sample rate(200)
% csi_filter = lowpass(input,Ws1,Wp1,Wp2,Ws2,fk)

%% subcarrier selection

[selected_h,selection] = subcarrier_selection_hb(Fs_h,csi_filter,10);

figure();   
subplot(211);
    plot(ts,selected_h);
    title('subcarrier selected heartbeat signal');
    xlabel('time/20ms');
    ylabel('amplitude');
subplot(212);
    plot(ts,csi_filter(:,selection));
    title('raw csi');
    xlabel('time/20ms');
    ylabel('amplitude');


%% filter

% csi_filter=filter(bandpass_hb,csi_raw);

figure();   
subplot(211);
    plot(ts,csi_filter);
    title('filtered heartbeat signal by bandpass');
    xlabel('time/20ms');
    ylabel('amplitude');
subplot(212);
    plot(ts,csi_raw);
    title('raw csi');
    xlabel('time/20ms');
    ylabel('amplitude');

%% EMD
% only For heartbeat singals, select the 2nd (and 3rd) IMF for further processing

emd_csi_h = emd_hb(selected_h,2);
%plot the IMFs
emd(selected_h(:,1),'Interpolation','pchip','display',0);

figure();
%     subplot(211);
%         plot(hamed_csi_h);
%         xlabel('time/20ms');
%         ylabel('amplitude');
%         title('raw heartbeat signal');
%     subplot(212);
        plot(ts,emd_csi_h);
        xlabel('time/20ms');
        ylabel('amplitude');
        title('filtered heartbeat signal by EMD');


%% Wavelet denoising

% heartbeat singal decomposition
% wavelet_csi_h = [];
% for i=1:num_subcarriers
%     [ch,lh]=wavedec(emd_csi_h(:,i),7,'sym7');
% % extract all 5 levels of approximate and detailed coeficients 
%     ah1=appcoef(ch,lh,'sym7',1);
%     dh1=detcoef(ch,lh,1);
%     ah2=appcoef(ch,lh,'sym7',2);
%     dh2=detcoef(ch,lh,2);
%     ah3=appcoef(ch,lh,'sym7',3);
%     dh3=detcoef(ch,lh,3);
%     ah4=appcoef(ch,lh,'sym7',4);
%     dh4=detcoef(ch,lh,4);
%     ah5=appcoef(ch,lh,'sym7',5);
%     dh5=detcoef(ch,lh,5);
%  % reconstruct signals based on the 5th detailed coefficients
%     dh1=zeros(size(dh1));
%     dh2=zeros(size(dh2));
%     dh3=zeros(size(dh3));
%     dh4=zeros(size(dh4));
%     ah5=zeros(size(ah5));
%     ch1=[ah5' dh5' dh4' dh3' dh2' dh1'];        %dh5 is non-zero
%     wavelet_csi_h=[wavelet_csi_h,waverec(ch1,lh,'sym7')];
% end

for i = 1:length(selection)
    wavelet_csi_h(:,i) = wavelet_breathe(emd_csi_h(:,i),7,'coif4');
end

figure();  
    subplot(211);
        plot(ts,emd_csi_h);    
        grid on
        xlabel('time/20ms');
        ylabel('Amplitude');
        title('emded signal');
    subplot(212);
        plot(ts,wavelet_csi_h); 
        grid on
        xlabel('time/100ms');
        ylabel('Amplitude');
        title('filtered heartbeat signal by WT');

%%  Variance
csiVar = VarianceAnalysis(wavelet_csi_h,1);
figure();
subplot(211);
    plot(ts,wavelet_csi_h); 
    grid on
    xlabel('time/100ms');
    ylabel('Amplitude');
    title('filtered heartbeat signal by WT');
subplot(212);
    plot(ts,csiVar);
    xlabel('Time')
    ylabel('Amplitude')
    title("After Variance selection")
    grid on

% %% Identify the best subcarrier (only one)
% % based on mean absolute deviation (MAD)
% 
% % For heartbeat signal:
% MAD_h =zeros(1,num_subcarriers);
% 
% % calculate the MAD matrix for each subcarrier
% for i=1:num_subcarriers
%     sum1 = 0;
%     for j=1:num_pkg
%         sum1 = abs(wavelet_csi_h(j,i)-mean(wavelet_csi_h(:,i))) + sum1;
%     MAD_h(i) = sum1 / num_pkg;
%     end
% end
% 
% % Sort the 'MAD_h' sequence by descend order, 
% % 'MAD_sort_h' saves the sorted result, and the 'index_h' saves the corresponding position of the data in MAD_sort_h
% [MAD_sort_h, index_h] = sort(MAD_h,'descend');    
% % 'max_h' returns the one subcarrier index corresponding to the fourth? first? largest average absolute deviation
% max_index_h = index_h(1);              
% 
% % figure(8);
% %     plot(yuh2(:,max2)); 
% %     xlabel('time/20ms');
% %     ylabel('Amplitude');
% %     title('selected subcarrier for heartbeat signal'); 
% 

%% PCA Principal Component Analysis (PCA) on raw data.
[coeff,score,latent] = pca(csiVar);
% Xcentered = score*coeff';
pca_csi_h = score(:,1);
figure();
plot(ts,pca_csi_h);
xlabel('Time')
ylabel('Amplitude')
title("After PCA")
grid on

% Smooth
avgsmooth=smooth(pca_csi_h,'rloess',50);
avgsmooth=avgsmooth./(max(avgsmooth)-min(avgsmooth));
figure();
plot(ts,avgsmooth);
xlabel('Time(s)')
ylabel('Amplitude')
title("After PCA")
grid on


%% Nelder Mead method
% %     % Fit a periodic signal to a sinusoidal signal by Nelder Mead method.
% % 
% %     options = optimset('Display','iter');   % display the related information of itering process.
%     x0  = [1,4,1];     % Initial point for Nelder Mead searching.
% 
% %     [x,fval_b,exitflag_b,output_b] = fminsearch(@optsine_b,x0,options);  % start searching for the minimun of the function 'optsine_b' and the corresponding variable values.
%     [x,fval_b,exitflag_b,output_b] = fminsearch(@optsine_b,x0);
%     for i =1:198       % return the fitted sine signal.
%         yub4(i) = x(1)*sin(x(2)*i+x(3));
%     end
%     figure(9);
%         plot(yub4,'r','linewidth',2);
%         hold on;
%         plot(yub3,'b','linewidth',1);
%         hold off;
%         legend('original singal','fitted signal');
%         xlabel('time/100ms');
%         ylabel('Amplitude');
%         title('fitted sine signal of respiratory signal'); 
% 
% %     [x,fval_h,exitflag_h,output_h] = fminsearch(@optsine_h,x0,options);
%     [x,fval_h,exitflag_h,output_h] = fminsearch(@optsine_h,x0);
%     fit_sine_h = zeros(length_scs,1);
%     for i =1:length_scs
%         yuh4(i) = x(1)*sin(x(2)*i+x(3));
%     end
%     figure(10);
%         plot(yuh4);
%         hold on;
%         plot(yuh3,'r');
%         hold off;
%         xlabel('time/20ms');
%         ylabel('Amplitude');
%         title('fitted sine signal of heartbeat signal'); 

%% Spectrum analysis
% based on DFT

% For heartbeat signal:
N = 2^(4+nextpow2(num_pkg));
% remove the first and last 200 points before DFT
% clear_h_t=clear_h(200:800)-mean(clear_h(200:800));     %为什么要-平均值？
% clear_h_t=wavelet_csi_h(100:1900,:)-mean(wavelet_csi_h(100:1900,:));     %为什么要-平均值？
% X=fft(clear_h_t,N);
fft_sub=fft(avgsmooth,N);
fftshift_sub=fftshift(fft_sub); 
f=(-N/2:N/2-1).*Fs_h./N;

figure();     
plot(f,(abs(fftshift_sub).^2));
%xlim([0,3]);
xlabel('frequncy / Hz');
ylabel('amplitude');
title('DFT of heartbeat signal pca');
grid on;    
        
[peak_h,index_h]=max((abs(fftshift_sub(1:N/2)).^2));
Heart_rate = abs(f(index_h));   
if (Heart_rate<0.9||Heart_rate>2)
    fprintf('DFT: %fHz......Warning: Your heart rate is abnormal, please try again to confirm!\n',Heart_rate);
else
    fprintf('DFT: Heart rate: %fHz. \n', Heart_rate);
end


%% Estimation based on the average gap among all peaks

[a_heartrate] = rate_estimation(ts,avgsmooth,1);
% For heart rate:
% [pmax_h,lmax_h]=findpeaks(avgsmooth,'minpeakdistance',20,'MinPeakProminence',0.1);   % find peaks and their locations 
% Heart_rate  = Fs_h / mean(diff(lmax_h)); 

[~,~,~,hight] = findpeaks(avgsmooth,'MinPeakWidth',10);
[pks_peak, locs_peak] = findpeaks(avgsmooth,ts,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
[pks_valley, locs_valley] = findpeaks(-1.*avgsmooth,ts,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
locs=[locs_peak(1:min([length(locs_peak),length(locs_valley)]));locs_valley(1:min([length(locs_peak),length(locs_valley)]))];
Heart_rate  = 1 / (2.*mean(diff(locs))); 

figure();
    plot(ts,avgsmooth,'b');
    hold on
    scatter(locs_peak,pks_peak,'r*');
    hold on
    scatter(locs_valley,-1.*pks_valley,'y*');
    xlabel('time / 20ms');
    ylabel('amplitude');
    title('Peaks of heartbeat signal');

if (Heart_rate<0.9||Heart_rate>2)
    fprintf('PEAK GAP: Heart rate: %fHz......Warning: Your heart rate is abnormal, please try again to confirm!\n',Heart_rate);
else
    fprintf('PEAK GAP: Heart rate: %fHz. \n', Heart_rate);
end


