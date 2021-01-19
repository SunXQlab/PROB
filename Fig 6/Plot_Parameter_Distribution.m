
%%%% 频数直方图与正态分布密度函数（拟合）
figure,
data=BETA{3};
[mu,sigma]=normfit(data);
C=[-1:1e-2:1];
[y,x]=hist(data,1000);
bar(x,y,'FaceColor','r','EdgeColor','w');box off
xlim([mu-3*sigma,mu+3*sigma])
a2=axes;
ezplot(@(x)normpdf(x,mu,sigma),[mu-3*sigma,mu+3*sigma])
set(a2,'box','off','yaxislocation','right','color','none')
title 'Distribution density function (fitted)' 

data=BETA{3};
alpha=0.01:0.01:0.75;
CI=zeros(length(alpha),2);
ind=1;
for i=alpha
   [MUHAT,SIGMAHAT,MUCI,SIGMACI] = normfit(data,i);  % NORMFIT Parameter estimates and confidence intervals for normal data.s
    CI(ind,:)=MUCI;
    ind=ind+1;
end
 CI  
   
   
   
   
%    
% fun=@(p,x) p(1)./x.*exp(-(((x)-p(2))/p(3)).^2/2);%对数正态分布密度函数
% figure,
% n=1000;
% [y x]=hist(data,n);
% %这里因为没有数据所以用上边两句模拟x,y数据
% 
% bar(x,y,1);hold on;%根据x,y数据画柱状图
% [maxy ind]=max(y);
% 
% p=nlinfit(x,y,fun,[maxy*x(ind),(x(ind)),1]);%拟合
% %p(1)~和幅度有关    p(2)~mu    p(3)~sigma 
% 
% yfit=fun(p,x); %计算拟合曲线
% plot(x,yfit,'r','linewidth',1);
% 
% figure, histfit(data,1000)