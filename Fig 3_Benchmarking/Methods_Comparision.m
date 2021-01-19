
%%%% Methods Comparesion; Comparing DPP-based dynamic network with static sample data-based network 

% load X_stage.mat
clear S1
ExpData=scExp_Net_Ctime(1:size(scExp_Net_Ctime,1)-1,:);
% ExpData=(ExpData-min(ExpData,[],2))./(max(ExpData,[],2)-min(ExpData,[],2));

%%%%%%% PCC-based correlation method
S1=corrcoef(ExpData');  
S1=abs(S1);  % score matrix
S1(logical(eye(size(S1,1))))=0;

% %%%%%%% MI-based correlation method
% clear S2
% S2=zeros(size(ExpData,1),size(ExpData,1));
% for i=1:size(ExpData,1)
%     for j=1:size(ExpData,1)
%         S2(i,j)=VectorMI(ExpData(i,:)',ExpData(j,:)');
%     end
% end
% 
% S2=abs(S2);
% S2=S2./max(max(S2));
% S2(logical(eye(size(S2,1))))=0;

% %%%%%%% CMI-based correlation method
% clear S3
% S3=zeros(size(ExpData,1),size(ExpData,1));
% for i=1:size(ExpData,1)
%     for j=setdiff([1:size(ExpData,1)],i)
%         S3(i,j)=kernelcmi(ExpData(i,:),ExpData(j,:),ExpData(setdiff([1:size(ExpData,1)],[i j]),:));
%     end
% end
% S3=abs(S3);
% S3=S3./max(max(S3));
% S3(logical(eye(size(S3,1))))=0;

% %%%%%%% Linear Regression 
%   clear S4
%   S4=zeros(size(ExpData,1),size(ExpData,1));
% for i=1:size(ExpData,1)
%     x=ExpData;
%     x(i,:)=[];
%     b = regress(ExpData(i,:)',[ones(1,size(x,2))',x']);
%     b(1)=[];
%     S4(i,setdiff([1:size(ExpData,1)],i)) = b;
% %     S4(i,i) =1;
% end
% 
% S4=abs(S4);
% S4=S4./max(max(S4));
% S4(logical(eye(size(S4,1))))=0;

%%%%%%% LASSO Regression
clear S2
for i=1:size(ExpData,1)
     y=ExpData(i,:);
     x=ExpData;
     x(i,:)=[];
    [A,FitInfo] = lasso_Sun(x',y,'Alpha',1,'CV',10,'Standardize',1,'Intercept_Flag',0);  % Find the coefficients of a regularized linear regression model using 10-fold cross-validation and the elastic net method with Alpha = 0.75. Use the largest Lambda value such that the mean squared error (MSE) is within one standard error of the minimum MSE.
    idxLambdaMinMSE = FitInfo.IndexMinMSE;  %Predict exam scores for the test data. 
    coef = A(:,idxLambdaMinMSE);
    S2(i,setdiff([1:size(ExpData,1)],i))=coef;

end
S2=abs(S2);
S2=S2./max(max(S2));
S2(logical(eye(size(S2,1))))=0;

C(logical(eye(size(C,1))))=0;

%%%%% GENIE3
clear S3
% VIM = GENIE3(NetGene');
VIM = GENIE3(ExpData');
S3=VIM';

%%% ARACNe
clear S4
S4=csvread('F:\Clinical Gene expression network Project\Reversion\Codes\Benchmarking_LPS_scRNAseq_data\Results\net_ARACNe.csv',1,1);
%%% CLR
clear S5
S5=csvread('F:\Clinical Gene expression network Project\Reversion\Codes\Benchmarking_LPS_scRNAseq_data\Results\net_CLR.csv',1,1);
%%% mrnet
clear S6
S6=csvread('F:\Clinical Gene expression network Project\Reversion\Codes\Benchmarking_LPS_scRNAseq_data\Results\net_mrnet.csv',1,1);


%%% SCODE
clear S7
S7=csvread('F:\Clinical Gene expression network Project\Reversion\Codes\Benchmarking_LPS_scRNAseq_data\Results\net_SCODE.csv',1,1);

%%% LEAP
clear S8
S8=csvread('F:\Clinical Gene expression network Project\Reversion\Codes\Benchmarking_LPS_scRNAseq_data\Results\net_LEAP.csv',1,1);


figure,
plot_roc( abs(S1), C~=0 ,1, 1) ;hold on % PCOR
plot_roc( abs(S2), C~=0 ,1, 1); hold on  % LASSO Regression
plot_roc( abs(S3), C~=0 ,1, 1) ;hold on  % GENIE3
plot_roc( abs(S4), C~=0 ,1, 1) ;hold on  % ARACNe
plot_roc( abs(S5), C~=0 ,1, 1); hold on  % CLR
plot_roc( abs(S6), C~=0 ,1, 1); hold on  % MRNET
plot_roc( abs(S7), C~=0 ,1, 1); hold on  % SCODE  % color: [0.47,0.67,0.19]
plot_roc( abs(S8), C~=0 ,1, 1); hold on  % LEAP  % color: [0.47,0.67,0.19]
plot_roc(S_PROB, C~=0 ,1, 1); % PROB  % color: [0.83,0.37,0.17]  
hold on; plot(0:0.01:1,0:0.01:1,  '-.','color',[0.50,0.54,0.53], 'LineWidth', 1);

L=legend('PCOR','LASSO','GENIE3','ARACNe','CLR','MRNET','SCODE','LEAP','PPOB');
set(L,'FontSize',10,'location','SouthEast')

SS=cell(9);
SS{1}=S1;SS{2}=S2;SS{3}=S3;SS{4}=S4;SS{5}=S5;SS{6}=S6;SS{7}=S7;SS{8}=S8;SS{9}=S_PROB;
for i=1:9
    labels=reshape(C~=0,size(C,1)*size(C,2),1);
    scores=reshape(SS{i},size(SS{i},1)*size(SS{i},2),1);
    AUC(i) =scoreAUC(labels, scores);
end

AUC

	

figure,
b=bar(1:9,AUC,'FaceColor','flat','BarWidth',0.5);
c=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;	0.3010 0.7450 0.9330;0 0.8 0.8;0.5 0.4 0.6;0.6350 0.0780 0.1840];	
c=[repmat([0.3010 0.7450 0.9330], 8,1);0 0.4470 0.7410];

for i=1:9
b.CData(i,:) = c(i,:);
end
for i = 1:9
text(i-0.4,AUC(i)+0.02,num2str(roundn(AUC(i),-3)));
end
% xlabel('GRN methods','FontSize',20);
ylabel('AUC','FontSize',20)
set(gca,'xtick',[1:9])
set(gca,'xticklabel',{'PCOR','LASSO','GENIE3','ARACNe','CLR','MRNET','SCODE','LEAP','PPOB'},'FontSize',10)
axis([0 10 0 0.75])
rotateticklabel(gca,45)
set(gca,'FontSize',15)

