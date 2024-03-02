%% EMD分解画图程序源码(公开版)
% 运行该代码前，请务必安装时频域分析工具箱TFA_Toolbox，获取地址：http://www.khscience.cn/docs/index.php/2020/04/09/1/
% 本程序2020.04定型，后续更新都在完整版代码中
%  原文链接 https://zhuanlan.zhihu.com/p/69630638
%  本代码地址：http://www.khscience.cn/docs/index.php/2020/04/09/11/
fs = 100;  %采样频率
t = 0:1/fs:10;
x = sin(20*pi*t);
y = 2*sin(1*pi*t);
sig = x+y;
PlotEMDandFFT(sig,fs);  %该函数的源码获取方式见"代码说明.docx"

% 关于函数： function imf = PlotEMDandFFT(y,Fs)
% 信号EMD分解与各IMF分量频谱对照图
% 输入：
% y为待分解信号
% Fs为采样频率
% 输出：
% imf为经emd分解后的各imf分量值，该值与直接调用emd函数返回的数值相同

%% 关于完整版代码：
% 如果需要封装好的画图函数（PlotEMDandFFT.m 和Fb_FFT.m）的源码，可在下述连接（完整版）获取。
% 源码中包括了店主最新代码，其中还包括：
% 整合版EMD函数：整合了G-Rilling工具箱和MATLAB自带工具箱的EMD分解方法，傻瓜式调用
% EMD分解图绘制函数：只绘制信号EMD分解图（不画频谱图）的函数，适合不需要频谱分析的场景
% 演示EMD画图工具函数调用方法的demos。
% 更为丰富、详细的注释。
% 完整版代码：https://item.taobao.com/item.htm?spm=a2oq0.12575281.0.0.15d11deb3eNsbA&ft=t&id=641792836460
%% 完整版代码重要更新：
% 20200821 解决了MATLAB版本识别，并自动选择工具箱版本
% 20200820 修复了调用MATLAB官方emd函数时imf分量未包含res的问题
% 20200805 加入了整合版EMD函数：整合了G-Rilling工具箱和MATLAB自带工具箱的EMD分解方法
% 20200804 加入了只绘制信号EMD分解图（不画频谱图）的函数
% 20200410 初始版本（公开版版本）