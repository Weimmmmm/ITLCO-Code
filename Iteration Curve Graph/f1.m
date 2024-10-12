clear all
clc

load('D:\研二\小论文\实验\收敛曲线\ITLCO.mat','every_bestf1');
load('D:\研二\小论文\实验\收敛曲线\TLCO.mat', 'every_bestf2');
load('D:\研二\小论文\实验\收敛曲线\MSMA.mat', 'every_bestf3');
load('D:\研二\小论文\实验\收敛曲线\IECO.mat', 'every_bestf4');
load('D:\研二\小论文\实验\收敛曲线\IAOA.mat', 'every_bestf5');
load('D:\研二\小论文\实验\收敛曲线\ESO.mat', 'every_bestf6');
% %------------ITLCO----------------------%
fname=16
d1=reshape(every_bestf1(fname,1:1976),8,[]);
D1=min(d1,[],1);
yyyy1=log(D1);
plot(yyyy1,'Color',[0.04 0.09 0.27],'lineWidth',1.5)
% plot(yyyyy1,'r','lineWidth',1)
set(gca,'XTick',0:40:1000)
%  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
  set(gca,'XLim',[0 200])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
hold on
%-----------------------TLCO-------------------
d2=reshape(every_bestf2(fname,1:1976),8,[]);
D2=min(d2,[],1);
yyyyy2=log(D2);
plot(yyyyy2,'Color',[0.1 0.8 0.9],'lineWidth',1.5)
% plot(yyyyy2,'g','lineWidth',1)
set(gca,'XTick',0:40:1000)
%  set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
 set(gca,'XLim',[0 200])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
 hold on



% %---------------------MSMA------------
d3=reshape(every_bestf3(fname,1:247),1,[]);
D3=min(d3,[],1);
yyyyy3=log(D3);
plot(yyyyy3,'Color',[0.42 0.35 0.80],'lineWidth',1.5)
set(gca,'XTick',0:40:1000)
 set(gca,'XLim',[0 200])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
% %--------------------IECO-------------------------
d4=reshape(every_bestf4(fname,1:1482),6,[]);
D4=min(d4,[],1);
yyyyy4=log(D4);
plot(yyyyy4,'Color',[0.00 0.79 0.34],'lineWidth',1.5)
set(gca,'XTick',0:40:1000)
 set(gca,'XLim',[0 200])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
% %--------------------IAOA-----------------------
% d5=reshape(every_bestf5(28,67:666),3,[]);
d5=reshape(every_bestf5(fname,1:1235),5,[]);
D5=min(d5,[],1);
yyyyy5=log(D5);
% yyyyy5(1)=[8.74107688677474];
plot(yyyyy5,'Color',[1 0 0],'lineWidth',1.5)
set(gca,'XTick',0:40:1000)
 set(gca,'XLim',[0 200])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on
%------------------------ESO---------------------
d6=reshape(every_bestf6(fname,1:741),3,[]);
D6=min(d6,[],1);
yyyyy6=log(D6);
plot(yyyyy6,'Color',[1 0.5 1],'lineWidth',1.5)
set(gca,'XTick',0:40:1000)
set(gca,'XLim',[0 200])
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'})
%set(gca,'XTicklabel',{'0','0.25','0.50','0.75','1.0'})
hold on


legend('ITLCO','TLCO','MSMA','IECO','IAOA','ESO');
xlabel('FEs(\times10^5)');
ylabel('Fitness Value(log)');
title('F16');
% grid on;
% print('F:\收敛曲线图\g1.tif','F:\收敛曲线图\d1','-dtif','-r300')

