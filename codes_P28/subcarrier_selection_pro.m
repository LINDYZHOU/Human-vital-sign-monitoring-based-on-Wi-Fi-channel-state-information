function [output,selection,location,snr] = subcarrier_selection_pro(input,type,amount)
%subcarrier_selection select better sensetive subcarriers
% input: subcarrier matrix
% type: 0-breath, 1-heartbeat
% amount: the amount of selected subcarriers

%initialize selection register
selection=zeros(1,amount);
snr=zeros(1,amount);
location=zeros(1,amount);
[num_pkg,num_subcarriers]=size(input);
fs=200;
N=2^(nextpow2(num_pkg)+4);
f = fs*(-(N/2):(N/2)-1)/N;
% delet DC
input=detrend(input);

for i=1:num_subcarriers
    subcarrier=input(:,i);
    fft_sub=fft(subcarrier,N);
    fftshift_sub=fftshift(fft_sub);
    % Spectrum symmetry.
    %     %test plot
    %     l=length(subcarrier);
    %     n=-l/2:l/2-1;
    %     f=n*200/l;
    %     plot(f,abs(fft_sub));
    %     axis([-1,1,-inf,inf]);
    %         f=Fs*(-N/2:(N/2-1))/N;
    %         plot(f,abs(fft_sub));
    %         axis([-1,1,-inf,inf]);

    switch type
        case 0  %breath
            center=abs(fftshift_sub(floor(N/2-0.6*N/fs):floor(N/2+0.6*N/fs),:));
            f_center = f(floor(N/2-0.6*N/fs):floor(N/2+0.6*N/fs));

    % center: (vector) the energy of  center part (low frequency) of the
    % spectrum symmetry.
    % N/2: from the 0Hz.
    % 0.6*N/fs: to the 0.6Hz of the whole spectrum.           

        case 1  %heartbeat
            center=abs(fftshift_sub(floor(N/2+1*N/fs):floor(N/2+2.0*N/fs),:));
            f_center = f(floor(N/2+1*N/fs):floor(N/2+2.0*N/fs));
    end
        % Spectrum symmetry.
        %     %test plot
        %     l=length(subcarrier);
        %     n=-l/2:l/2-1;
        %     f=n*200/l;
        %     plot(f,abs(fft_sub));
        %     axis([-1,1,-inf,inf]);

    % Calculate SNR
    avg = mean(center);
    centerp=center.^2;
    SNR = centerp./avg^2;
    snr_max=max(SNR);
%     figure();     
%     plot(f_center,SNR);
% %     % plot(f,psd);
%     xlim([1,2]);
%     xlabel('frequncy / Hz');
%     ylabel('SNR');
% %     title('DFT of heartbeat signal pca');
%     grid on;  

    % Calculate peaks and those locations
    [pks,locs]=findpeaks(SNR);
    pk=sortrows([pks,locs],1,'descend');
    Hz0=(pk(1,2)+pk(2,2))/2;
    step=fs/N;
    pk(:,2)=(pk(:,2)-Hz0).*step;
    max_pk=max(pks);
    pks(pks>=max_pk)=0;
    max2_pk=max(pks);
    % plot(centerp);

    % Judge whether there is only one peak or not.
    % We only select those subcarriers with one peak to reach less interference.
    if max_pk > max2_pk*2
        % Compare the maxsnr of this subcarrier with existed mini snr
        [snr_min,loc_min]=min(snr);
        if snr_max > snr_min
            snr(loc_min)=snr_max;
            selection(loc_min)=i;
            location(loc_min)=abs(pk(1,2));
        end
    end
end

% Remove 0
selection(find(selection==0))=[];
location(find(location==0))=[];
snr(find(snr==0))=[];
output=input(:,selection);

%             figure();
%             plot(f_center,SNR);
%             title(['No.',num2str(i),' subcarrier']);

end

