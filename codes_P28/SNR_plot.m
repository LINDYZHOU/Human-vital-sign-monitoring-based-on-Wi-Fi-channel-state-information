%SNR plot express

% load csi data
[Rx1,Rx2,sequence]=CSIdatapro(ax20_regular60s_zxl,'hesu');

% frame to reduce the inputting data scale
%every 10 seconds is a frame
frame_length=6;
%every 10 seconds produce a new frame
frame_interval=6;
%sampling rate is 200Hz
rate=200;
% divide the data flow into frame and select the credible frame
Rx1_amp_framed=framing(Rx1{1},sequence,frame_length,frame_interval);
Rx2_amp_framed=framing(Rx2{1},sequence,frame_length,frame_interval);
% Rx1_pha_framed=framing(Rx1{2},sequence,frame_length,frame_interval);
% Rx2_pha_framed=framing(Rx2{2},sequence,frame_length,frame_interval);
% utilize phase difference instead of amplitude or phase to reduce the
% influnce from thermal noise
% phase_diff=Rx1{2}-Rx2{2};
% pha_diff_framed=framing(phase_diff,sequence,frame_length,frame_interval);

% calculate SNR of each subcarrier every frame
for input=1:2
    switch input
        case 1
            csi_framed=Rx1_amp_framed;
        case 2
            csi_framed=Rx2_amp_framed;
    end
    num_frame = length(csi_framed);
    for n = 1:num_frame
        frame = csi_framed{n};
        [num_pkg,num_subcarriers]=size(frame);

        inter_power=abs(mean(frame-hampel(frame),1));
        frame=hampel(frame);
        s_power=mean(frame,1);
        diff=frame-ones(num_pkg,1)*s_power;
        d_power=max(diff,[],1);

        ssnr=d_power.^2./(s_power+inter_power).^2;
        ssnr_frame(n,:)=ssnr;


    end
    var_ssnr=var(ssnr_frame);
    var_snr_avg=mean(var_ssnr);
    [~,index]=sortrows(var_ssnr','descend');

    figure();
    num_sub=5;
    for num=1:num_sub
    plot(1:6:num_frame*6,ssnr_frame(:,index(num)),'DisplayName',['第',num2str(index(num)),'个子载波']);
    hold on
    end
    hold off
    title(['Rx',num2str(input),' amplitude']);
    xlabel('时间/s');
    ylabel('SSNR');


end













