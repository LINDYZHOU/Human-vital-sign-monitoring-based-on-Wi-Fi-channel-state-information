function [output,selection] = subcarrier_selection_h(fs,input,threshold)
%subcarrier_selection_h select better sensetive subcarriers for heartrate
%estimation.
% Zhou Xiaolin 20220426
% fs: sampling rate of input data
% input: subcarrier matrix
% threshold: the ratio between the peaks of different frequency components
% and the average.

[num_pkg,num_subcarriers] = size(input);
selection=[];
output=[];
% preparation for FFT.
N = 2^(4+nextpow2(num_pkg));
f = fs*(-(N/2):(N/2)-1)/N;

finded = 0;

while(finded == 0)
    
    for i=1:num_subcarriers
        subcarrier=input(:,i);
        subcarrier=detrend(subcarrier); % delet DC
        fft_sub=fft(subcarrier,N);
        fftshift_sub=fftshift(fft_sub); 
        % Spectrum symmetry.
        %test plot
%         f=N*20/l;
%         N = 2^nextpow2(num_pkg);
%         f = Fs*(0:(N/2))/N;
%         plot(f,abs(fftshift_sub));
%         axis([-1,1,-inf,inf]);

        center=abs(fftshift_sub(floor(N/2+1*N/fs):floor(N/2+2.0*N/fs),:));
        f_center = f(floor(N/2+1*N/fs):floor(N/2+2.0*N/fs));
        % center: (vector) the energy of  center part (low frequency) of the
        % spectrum symmetry.
        % N/2+1*N/fs: from the 1Hz.
        % N/2+2.0*N/fs: to the 2Hz.
    
        avg = mean(center);
        centerp = center.^2;
        SNR = (center.^2)./avg^2;
    
        if sum(SNR>threshold) ~= 0
    %     if sum(center>avg*threshold) ~= 0
            % when there is SNR bigger than threshold in
            % vector center, this subcarrier has good performance in low
            % frequency, thus selected.

            [pks,locs]=findpeaks(centerp);
            pk=sortrows([pks,locs],1,'descend');
            Hz0=(pk(1,2)+pk(2,2))/2;
            step=2/(length(centerp)-1);
            pk(:,2)=(pk(:,2)-Hz0).*step;
            max_pk=max(pks);
            pks(pks>=max_pk)=0;
            max2_pk=max(pks);
            % plot(centerp);
            % only select those subcarriers with one-peak SNR.

            if max_pk > max2_pk*2
                selection=[selection,i];
                output=[output,subcarrier];
%                 location=[location,abs(pk(1,2))];

%             figure();
%             plot(f_center,SNR);
%             title(['No.',num2str(i),' subcarrier']);

            end

        end
    end

    if isempty(selection)
%         error("no subcarrier is suitable!");
    threshold = (0.95)*threshold;
    else
         finded = 1;
    end

    if threshold < 1
        error("no subcarrier is suitable!");
    end

end
end

