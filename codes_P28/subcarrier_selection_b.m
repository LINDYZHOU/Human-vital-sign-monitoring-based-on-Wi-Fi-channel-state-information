function [output,selection,location,threshold] = subcarrier_selection_b(input,threshold)
%subcarrier_selection select better sensetive subcarriers
% input: subcarrier matrix
% threshold: the ratio between the peaks of different frequency components
% and the average.
flag=1;
while(flag)
    [num_pkg,num_subcarriers]=size(input);
    selection=[];
    output=[];
    location=[];
    SNR_max = [];
    fs=200;
    N=2^(nextpow2(num_pkg)+2);

    for i=1:num_subcarriers
        subcarrier=input(:,i);
%         subcarrier=detrend(subcarrier); % delet DC
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
        center=abs(fftshift_sub(floor(N/2-N/fs): ...
            floor(N/2+N/fs),:));
        % center: (vector) the energy of  center part (low frequency) of the
        % spectrum symmetry.
        % num_pkg/2: from the 0.
        % num_pkg/20: to the 1/20 of the whole spectrum.

        avg = mean(center);
        centerp = center.^2;
        SNR = centerp./(avg^2);

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
            if max_pk > max2_pk*2
                selection=[selection,i];
                output=[output,subcarrier];
                location=[location,abs(pk(1,2))];
                SNR_max = [SNR_max,max(SNR)];

%             figure();
%             plot(f_center,SNR);
%             title(['No.',num2str(i),' subcarrier']);

            end
        end
    end

    if length(selection) < 3
        %error("no subcarrier is suitable!");
        threshold=threshold*0.95;
    else
        flag=0;
    end

    if length(selection) > 10
    % too many selected SCs
    [SNR_sort_max, index_h] = sort(SNR_max,'descend');    
    % Sort the 'SNR_max' sequence by descend order, 
    % 'SNR_sort_max' saves the sorted result
    % the 'index_h' saves the corresponding position of the data in SNR_sort_max
    SNR_sort_index = index_h(1:10);
    % 'SNR_sort_index' returns the FIRST 10 subcarrier index corresponding to 
    % the  largest average SNR
    selection=selection(SNR_sort_index);
    output=output(:,SNR_sort_index);
    location=location(SNR_sort_index);
    end
    
end
end

