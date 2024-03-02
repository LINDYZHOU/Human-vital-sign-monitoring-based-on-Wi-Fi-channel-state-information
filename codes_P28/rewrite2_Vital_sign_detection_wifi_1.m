%% Notes
%   Update on April 16, 2021 by Wang Linguo
%   The function of this program is to achieve the estimation of respiratory rate, extraction of heartbeat signal from respiratory signal and the estmation  of heart rate.
%   Fs = 50Hz, Sampling time t = 20s

        
%% Read the CSI data from .dat file

% length_scs is the data length of each subcarrier, which should be 1000. 
% Considering there might be uncertain unmber of packges lost, set 'lx' = 990.
num_pkg = 990;
Fs_h=50;
[rx1_csi,rx2_csi,rx3_csi] = readCSIn ('E:\FYP\220307_Wang Linguo-Vital Sign Detection\data\20210304\sit100_LOS\221103033220210304.dat',num_pkg);

% 'data' is a matrix of 990*90 including all 90 subcarriers data
data = [rx1_csi;rx2_csi;rx3_csi]';  
[num_pkg,num_subcarriers]=size(data);


%% Hampel filter  
% Remove outliners and replace them with local median    

% filter heartbeat signal    
[hamed_csi_h, outlier_h]=hampel(data,3);     %the length of window of both side is 3

figure();   
    subplot(211);
        plot(hamed_csi_h);
        title('filtered heartbeat signal by Hampel');
        xlabel('time/20ms');
        ylabel('amplitude');
    subplot(212);
        plot(selected_h);
        title('heartbeat signal selected by SNR');
        xlabel('time/20ms');
        ylabel('amplitude');


%% subcarrier selection

[selected_h,selection] = subcarrier_selection(Fs_h,hamed_csi_h,10);

%% EMD
% only For heartbeat singals, select the 3rd and 4th IMF for further processing

% decompose heartbeat signal
emd_csi_h = zeros(num_pkg,length(selection));
for i = 1:length(selection)        %each subcarrier
    [imf_h,residual_h] = emd(selected_h(:,i),'Interpolation','pchip','display',0);
    %Name-'interpolation', consisting of  either 'spline' or 'pchip'.
    % 'spline', if x is a smooth signal
    % 'pchip', if x is a nonsmooth signal
    emd_csi_h(:,i)=imf_h(:,4)+imf_h(:,3);    

    %frequency domain
%     figure();
%     N = 2^nextpow2(num_pkg);
%     f = Fs_h*(-(N/2):(N/2)-1)/N;
%     imf_dc=detrend(imf_h(:,4)); % delet DC
%     fft_imf=fft(imf_dc,N);
%     fftshift_imf=fftshift(fft_imf); 
%     plot(f,abs(fftshift_imf));

end

%plot the IMFs
emd(selected_h(:,3),'Interpolation','pchip','display',0);

%piot the IMFs and FFT
PlotEMDandFFT(selected_h(:,2),Fs_h);

figure();
%     subplot(211);
%         plot(hamed_csi_h);
%         xlabel('time/20ms');
%         ylabel('amplitude');
%         title('raw heartbeat signal');
%     subplot(212);
        plot(emd_csi_h(:,3));
        xlabel('time/20ms');
        ylabel('amplitude');
        title('filtered heartbeat signal by EMD');


%% Wavelet denoising

% heartbeat singal decomposition
wavelet_csi_h = zeros(num_pkg,num_subcarriers);
for i=1:num_subcarriers
    [ch,lh]=wavedec(emd_csi_h(:,i),5,'sym7');

% extract all 5 levels of approximate and detailed coeficients 
    ah1=appcoef(ch,lh,'sym7',1);
    dh1=detcoef(ch,lh,1);
    ah2=appcoef(ch,lh,'sym7',2);
    dh2=detcoef(ch,lh,2);
    ah3=appcoef(ch,lh,'sym7',3);
    dh3=detcoef(ch,lh,3);
    ah4=appcoef(ch,lh,'sym7',4);
    dh4=detcoef(ch,lh,4);
    ah5=appcoef(ch,lh,'sym7',5);
    dh5=detcoef(ch,lh,5);

 % reconstruct signals based on the 5th detailed coefficients
    dh1=zeros(size(dh1));
    dh2=zeros(size(dh2));
    dh3=zeros(size(dh3));
    dh4=zeros(size(dh4));
    ah5=zeros(size(ah5));
    ch1=[ah5' dh5' dh4' dh3' dh2' dh1'];        %dh5 is non-zero
    wavelet_csi_h(:,i)=waverec(ch1,lh,'sym7');
end
% figure(6);  
%     subplot(211);
%         plot(yuh1);           
%         xlabel('time/20ms');
%         ylabel('Amplitude');
%         title('raw signal');
%     subplot(212);
%         plot(yuh2); 
%         xlabel('time/100ms');
%         ylabel('Amplitude');
%         title('filtered heartbeat signal by WT');


%% Identify the best subcarrier (only one)
% based on mean absolute deviation (MAD)

% For heartbeat signal:
MAD_h =zeros(1,num_subcarriers);

% calculate the MAD matrix for each subcarrier
for i=1:num_subcarriers
    sum1 = 0;
    for j=1:num_pkg
        sum1 = abs(wavelet_csi_h(j,i)-mean(wavelet_csi_h(:,i))) + sum1;
    MAD_h(i) = sum1 / num_pkg;
    end
end

% Sort the 'MAD_h' sequence by descend order, 
% 'MAD_sort_h' saves the sorted result, and the 'index_h' saves the corresponding position of the data in MAD_sort_h
[MAD_sort_h, index_h] = sort(MAD_h,'descend');    
% 'max_h' returns the one subcarrier index corresponding to the fourth? first? largest average absolute deviation
max_index_h = index_h(1);              

% figure(8);
%     plot(yuh2(:,max2)); 
%     xlabel('time/20ms');
%     ylabel('Amplitude');
%     title('selected subcarrier for heartbeat signal'); 


%%  Normolize the signal and remove the DC component 

max_value = max(abs(wavelet_csi_h(:,max_index_h)));      % find the max value of the MAD subcarrier

% normolization
for i=1:num_pkg
   normolized_csi_h(i,max_index_h) = wavelet_csi_h(i,max_index_h) / max_value;
end

global clear_h_g;       % used in function: optsine_h
clear_h_g=normolized_csi_h(:,max_index_h)-mean(normolized_csi_h(:,max_index_h));       % remove the DC component
clear_h = clear_h_g;


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
N = 2^nextpow2(num_pkg);
% remove the first and last 200 points before DFT
% clear_h_t=clear_h(200:800)-mean(clear_h(200:800));     %为什么要-平均值？
clear_h_t=emd_csi_h(100:900,3)-mean(emd_csi_h(100:900,3));     %为什么要-平均值？
% X=fft(clear_h_t,N);
X=fft(clear_h_t,N);
n=0:N/2;
f=n.*Fs_h./N;

figure();     
plot(f,(abs(X(1:N/2+1)).^2));
xlim([0,3]);
xlabel('frequncy / Hz');
ylabel('amplitude');
title('DFT of heartbeat signal 3');
grid on;    
        
[peak_h,index_h]=max((abs(X(1:N/2+1)).^2));
Heart_rate = f(index_h);   
if (Heart_rate<0.9||Heart_rate>2)
    fprintf('DFT: %fHz......Warning: Your heart rate is abnormal, please try again to confirm!\n',Heart_rate);
else
    fprintf('DFT: Heart rate: %fHz. \n', Heart_rate);
end


%% Estimation based on the average gap among all peaks

% For heart rate:
[pmax_h,lmax_h]=findpeaks(clear_h_g,'minpeakdistance',20,'MinPeakProminence',0.1);   % find peaks and their locations 
Heart_rate  = Fs_h / mean(diff(lmax_h)); 
% figure(14);
%     plot(yuh3,'b');
%     hold on;
%     scatter(lmax_h,pmax_h,'r*');
%     hold off;
%     xlabel('time / 20ms');
%     ylabel('amplitude');
%     title('Peaks of heartbeat signal');
if (Heart_rate<0.9||Heart_rate>2)
    fprintf('PEAK GAP: Heart rate: %fHz......Warning: Your heart rate is abnormal, please try again to confirm!\n',Heart_rate);
else
    fprintf('PEAK GAP: Heart rate: %fHz. \n', Heart_rate);
end


%% subfunctions

% define the function which is used to measure the  minimum mean square error between 
% fitted sine signal and the input respiratory signal.
% function output = optsine_b(ab)  
%     output = 0;
%     global yub3;
%     for i=1:198
%         sum2 =  (ab(1)*sin(ab(2)*i+ab(3))-yub3(i))^2;
%         output = output + sum2;
%     end
% end

% define the function which is used to measure the  minimum mean square error between fitted sine signal and the input heartbeat signal.
function output = optsine_h(ab)
    output = 0;
    global yuh3;
    for i=1:length_scs
        sum2 =  (ab(1)*sin(ab(2)*i+ab(3))-yuh3(i))^2;
        output = output + sum2;
    end
end