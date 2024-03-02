function [output] = pp_estimation(ts,input,type)
%this function is used to do respiratory/heartbeat peak-valley rate estimation.
%input: the smoothed data 
%type: "0"-breathe, "1"-heartbeat
%output: the estimated rate.
frame_length = 10;
[~,~,~,hight] = findpeaks(input,'MinPeakWidth',10);
[pks_peak, locs_peak] = findpeaks(input,ts,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
[pks_valley, locs_valley] = findpeaks(-1.*input,ts,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
switch type
    case 0      %breath rate
        %if there is 2 breaths in 1 frame
        if length(pks_peak) >= 2 && length(pks_valley) >= 2
            % no less than 2 times
            %combine location of peaks and valleys
            locs=[locs_peak(1:min([length(locs_peak),length(locs_valley)]));locs_valley(1:min([length(locs_peak),length(locs_valley)]))];
            %find the difference time of the adjacent peak and valley, which is
            %half of the duration
            Interval = abs(diff(locs,1,1));
            duration = mean(Interval);
            output = 1/(2*duration);
        else
            % less than 2 times
            % mirror the avgsmmoth twice  
            wave=[fliplr(input')';input;fliplr(input')'];
        %     figure();
        %     plot(breath_wave);
            % 3 times longer and then find peaks
            % make man-made duration to hlep findpeaks
            [~,~,~,hight] = findpeaks(wave,'MinPeakWidth',10);
            pks = findpeaks(wave,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
            vals = findpeaks(-1.*wave,'MinPeakProminence',0.2*median(hight),'Annotate','extents');
            breath_times=(length(pks)+length(vals))/6;
            output=breath_times/frame_length;
    end
    % text(6,0,['respiratory rate:',num2str(respiratory_rate(n))])
    %str = {"breathtimes = ",num2str(breathtimes(i)),"respiratory interval: (seconds)",num2str(peakInterval),"respiratory duration = (seconds)", num2str(respiratory_duration),"respiratory rate = (Hz)", num2str(respiratory_rate(i))};
    case 1      %heart rate
        locs=[locs_peak(1:min([length(locs_peak),length(locs_valley)]));locs_valley(1:min([length(locs_peak),length(locs_valley)]))];
        Interval = abs(diff(locs,1,1));
        duration = mean(Interval);
        output  = 1 / (2*duration); 
end
% figure();
% % findpeaks(avgsmooth,t_sampling,'MinPeakProminence',2,'Annotate','extents');
% plot(t_sampling,avgsmooth);
% % text(1,1,str)
% % grid on
% xlabel('Time')
% ylabel('Amplitude')
% title(["Frame ";num2str(n)])
% figure();
%     plot(ts,input,'b');
%     hold on
%     scatter(locs_peak,pks_peak,'r*');
%     hold on
%     scatter(locs_valley,-1.*pks_valley,'y*');
%     xlabel('时间');
%     ylabel('幅度');
%     title('重构心跳信号波峰、波谷示意图');
end