%% B&A subcarriers selection    
figure();
    subplot(311)
    plot(t,frame);
    grid on
    xlabel('时间')
    ylabel('幅度')
    title("原始数据")
    subplot(312);
    plot(t,csi_selected_h); 
    xlabel('时间')
    ylabel('幅度')
    title('SSNR子载波选择后');
    grid on
    subplot(313);
    plot(t,csi_var_selected_h); 
    xlabel('时间')
    ylabel('幅度')
    title('最大方差子载波选择后');
    grid on
%% EMD, coif4
figure();
    subplot(211)
    plot(t,emd_csi_h);
    grid on
    xlabel('时间')
    ylabel('幅度')
    title("EMD后")
    subplot(212);
    plot(t,csi_wavelet_h); 
    xlabel('时间')
    ylabel('幅度')
    title('EMD+coif4 DWT后');
    grid on
    %% db4
    figure();
    subplot(211)
    plot(t_db4,csi_wavelet_h_db4);
    grid on
    xlabel('时间')
    ylabel('幅度')
    title("db4 DWT后")
    subplot(212);
    plot(t,csi_wavelet_h); 
    xlabel('时间')
    ylabel('幅度')
    title('EMD+coif4 DWT后');
    grid on
