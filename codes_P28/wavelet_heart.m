function output = wavelet_heart(input,order,wvname)
% Heartbeat Wavelet Transform
% 用coif4小波对原始信号进行7层分解并提取系数
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
% d_l = wrcoef('d',c,l,wvname,order-1);
% a_l = wrcoef('a',c,l,wvname,order-1);
%WRCOEF 重构近似/细节系数
%X = WRCOEF('type',C,L,'wname',N)
%Argument 'type' determines whether approximation ('type' = 'a') or 
% detail ('type' = 'd') coefficients are reconstructed.
output = d_h ;
% [d_l,d_h] = detcoef(c,l,[order-1,order]);
%  DETCOEF提取一维细节系数。
%  D = DETCOEF(C,L,N)从小波分解结构[C,L]中提取出N级的细节系数(the detail coefficients)。
%  Level N must be an integer such that 1 <= N <= NMAX where NMAX = length(L)-2.
%  D = DETCOEF(C,L)提取最后一级NMAX的细节系数。
%  If N is a vector of integers such that 1 <= N(j) <= NMAX: 