%画柱状图

Y=[89.93,94.74                 %2组数据
    89.56,90.07
    ];
figure();
set(gcf,'position',[100,100,800,500])   %调整图的位置
X=1:size(Y,1);
h=bar(X,Y,0.6);       %画出两个直方图，宽度0.9，可调
set(gca,'XTickLabel',{'IEEE 802.11ac','IEEE 802.11ax'},'FontSize',16,'FontName','Times New Roman'); 
%修改横坐标名称、字体
set(h(1),'FaceColor',[1 0.4 0])     % 设置条形图颜色，图1
set(h(2),'FaceColor',[0 0.7 1])     % 设置条形图颜色，图2
ylim([0,105]);      %y轴刻度
%修改x,y轴标签，中英文字体分开
ylabel('\fontname{宋体}准确率/%');
xlabel('\fontname{宋体}不同协议'); 
%修改图例，中英文字体分开
legend({'\fontname{宋体}坐','\fontname{宋体}卧'}, 'FontSize',14,'Location','northeastoutside');
set(gca,'xtick',1:4);   %x轴刻度
Y_1=roundn(Y,-2);   %调整y轴数字的精度，保留小数点后几位
%在柱状图上标数字（百度找的，哈哈哈，出处忘了，sorry），距离可调
for i = 1:length(X)
    text(X(i)-0.15,Y_1(i,1),num2str(Y_1(i,1)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14,'FontName','Times New Roman');
    text(X(i)+0.15,Y_1(i,2),num2str(Y_1(i,2)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14,'FontName','Times New Roman');
end