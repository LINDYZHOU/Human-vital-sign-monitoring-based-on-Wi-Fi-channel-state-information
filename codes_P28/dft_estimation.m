function [output] = dft_estimation(input,rate)
% This function is for DFT rate estimation.
%input: the processed signal
%rate: 20-respiratory, 200-heartrate
%output: estimated rate 

    N = 2^(4+nextpow2(length(input)));
    fft_sub=fft(input,N);
    fftshift_sub=fftshift(fft_sub); 
    f=(-N/2:N/2-1).*rate./N;
    
%     figure();     
%     plot(f,(abs(fftshift_sub).^2));
%     % plot(f,psd);
%     xlim([0,1]);
%     xlabel('frequncy / Hz');
%     ylabel('power');
%     title('DFT of respiration signal ');
%     grid on;  
    
    psd = (abs(fftshift_sub(1:N)).^2)./N;
    [~,index_h]=max(psd);
    output =abs(f(index_h));   
    
end