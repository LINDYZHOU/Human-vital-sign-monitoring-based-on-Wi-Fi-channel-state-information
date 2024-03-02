function output = wavelet_breathe(input,order,wvname,frequency)
% Fast/Slow Breath Wavelet Transform
%对原始信号进行小波分解
% 用db6小波对原始信号进行6层分解并提取系数
[c,l]=wavedec(input,order,wvname);

% function [c,l] = wavedec(x,n,IN3,IN4)：多层一维小波分解。
% WAVEDEC 使用一个特定的小波“wname”或一组特定的小波分解滤波器进行多层一维小波分析。
% [C,L] = WAVEDEC(X,N，'wname') 使用'wname'返回信号X在N级的小波分解。
% 输出分解结构包含小波 the wavelet decomposition vector C和 the bookkeeping vector L
% N 必须是一个严格的正整数。

% 对[c,l]提取近似系数
% a_coef=appcoef(c,l,wvname,order);
% APPCOEF: 提取一维小波变换近似系数。
% A = APPCOEF(C,L,'wname',N) 使用小波分解结构计算N级的近似系数[C,L]
% Level N must be an integer such that 0 <= N <= length(L)-2.
% A = APPCOEF(C,L,'wname') extracts（提取） the approximation coefficients（近似系数）
% at the last level length(L)-2.

% 对[c,l]提取细节系数
d_h = wrcoef('d',c,l,wvname,order);
d_l=wrcoef('d',c,l,wvname,order-1);
a_h = wrcoef('a',c,l,wvname,order);
%WRCOEF 重构近似/细节系数
%X = WRCOEF('type',C,L,'wname',N)
%Argument 'type' determines whether approximation ('type' = 'a') or
% detail ('type' = 'd') coefficients are reconstructed.
if abs(frequency) > 0.3
    output = d_l;
else
    output = d_h;
end

% output = a_h;

end

% [d_l,d_h] = detcoef(c,l,[order-1,order]);
%  DETCOEF提取一维细节系数。
%  D = DETCOEF(C,L,N)从小波分解结构[C,L]中提取出N级的细节系数(the detail coefficients)。
%  Level N must be an integer such that 1 <= N <= NMAX where NMAX = length(L)-2.
%  D = DETCOEF(C,L)提取最后一级NMAX的细节系数。
%  If N is a vector of integers such that 1 <= N(j) <= NMAX:


% 强制消噪
%对信号进行强制性消噪处理并图示结果
% c1=zeros(1,length(c)-length(d_h)-length(d_l));
% 将需要的细节系数与其他系数（均为0）合并
% c1=[c1 d_l d_h]; % c1有用信息仅来自细节系数d_l d_h
% output=waverec(c1,l,wvname);% 用c1来重构c
% WAVEREC：多层一维小波重建。与WAVEDEC相对应
% WAVEREC使用特定的小波('wname'，参见WFILTERS)或特定的重建滤波器(Lo_R和Hi_R)来执行多层一维小波重构。
% X = WAVEREC(C,L，'wname')根据多级小波分解结构对信号X进行重构


% % 默认阈值对信号进行消噪处理并图示结果
% % 1、用ddencmp函数获得信号的默认阈值
% [thr,sorh,keepapp]=ddencmp('den','wv',s);
% %使用ddencmp函数的去噪功能对原始信号s进行小波分解，并得到阈值thr，
% % 阈值方式，是否允许近似保持系数
%
% %  DDENCMP: 去噪或压缩的默认值
% %  [THR,SORH,KEEPAPP,CRIT] = DDENCMP(IN1,IN2,X)：
% % 返回输入向量或矩阵X的默认去噪或压缩值，使用小波或小波包，可以是一维或二维信号。其中：
% %  THR是阈值，SORH是函数（自动）选择的阈值方式：
% % 软阈值（SORH = s）或硬阈值（SORH = h），KEEPAPP允许保持近似系数（0、1，
% % CRIT(仅用于小波包)是熵名(参见WENTROPY)。
% %  IN1是den（去噪）或者cmp（压缩）
% %  IN2是wv（小波分解）或者wp（小波包分解）。
%
% s2=wdencmp('gbl',c,l,'db1',3,thr,sorh,keepapp);
% %利用ddencmp得到的阈值、阈值方式、近似保持系数对小波分解后的信号c进行去噪
%
% % WDENCMP使用小波对信号或图像进行去噪或压缩处理。
% % [XC,CXC,LXC,PERF0,PERFL2] = WDENCMP('gbl'，X，'wname'，N,THR,SORH,KEEPAPP)返回输入信号X 的去噪或压缩版本的XC，使用全局正阈值THR进行小波系数阈值化。
% %  附加输出参数[CXC,LXC]是XC的小波分解结构，
% % PERFL2和PERF0是以百分比表示的L^2恢复和压缩分数。
% % PERFL2 = 100*(CXC的向量范数/ C的向量范数)^2其中[C,L]为X的小波分解结构。
% % 小波分解在第N层执行，'wname'是一个包含小波名称的字符串。
% % SORH ('s'或'h')用于软阈值或硬阈值(有关详细信息，请参阅WTHRESH)。
% % 如果KEEPAPP = 1，逼近系数不能设置阈值，否则是可能的。
%
%
% subplot(2,2,3);
% plot(s2);
% title('默认阈值消噪后的信号');
% xlabel('样本序号n');
% ylabel('幅值A');
