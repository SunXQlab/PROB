%%%% Fig 2-associated code

scExpression=csvread('F:\Clinical Gene expression network Project\Reversion\Data\LPS_WT_scExpression.csv',1,1);
Capture_time=csvread('F:\Clinical Gene expression network Project\Reversion\Data\LPS_Capture_time.csv',1,1);

TF_network_ind=csvread('F:\Clinical Gene expression network Project\Reversion\Data\TF_network_ind.csv',1,1);
TF_network=csvread('F:\Clinical Gene expression network Project\Reversion\Data\TF_network.csv',1,1);

input_data=[scExpression;Capture_time'];
cd('F:\Clinical Gene expression network Project\Reversion\Codes\Benchmarking_LPS_scRNAseq_data')
save scExpression_Time.mat 
[val,ind]=sort(Capture_time);
input_data=input_data(:,ind);
save input_data.mat input_data
tic
[Data_smooth,DPT,DPP]=Progression_Inferrence(input_data);
toc
NetGene=Data_smooth(TF_network_ind,:);
% [val,ind]=sort(DPT);
% output_data=input_data(:,ind);
% NetGene=output_data(TF_network_ind,:);
% DPT_sort=val';
tic
[Para_Post_pdf,S]=ODE_BayesianLasso(NetGene,DPP);  %DPT_sort
toc
   save DPT.mat DPT
   save DPP.mat DPP
%    scExp_Net_Ctime=[NetGene;val'];
scExp_Net_Ctime=[NetGene;DPP];
   save scExp_Net_Ctime.mat scExp_Net_Ctime
   
   C=TF_network;
   clear S_PROB
   for i=1:size(C,1)
        S_PROB(i,i)=0;
        S_PROB(i,setdiff(1:size(C,1),i))=S(i,1:end-1);
   end
  save S_PROB.mat  S_PROB
  
  %%% Real time vs DPT
  % Spearman correlation
[RHO,PVAL]=corr(input_data(size(input_data,1),:)',DPT,'type','Spearman')
figure,
scatter(input_data(size(input_data,1),:)/max(input_data(size(input_data,1),:)),DPT/max(DPT),40,DPT/max(DPT),'fill')
set(gca,'FontSize',15)
xlabel('Real progression','FontSize',20);
ylabel('Inferred pseudo-progression','FontSize',20)
colorbar
hold on, plot(0:0.01:1,0:0.01:1,'-.','Color',[0.5,0.5,0.5])
title(['rho=',num2str(RHO)])

% R square
mdl = fitlm(input_data(size(input_data,1),:)',DPT);
mdl.Rsquared.Ordinary  % R^2

%% For cytoscape visualization
clear Mean_Para
Act_Inh_strength=zeros(size(S));Act_Inh_strength=zeros(size(S));
for i=1:size(S,1)
Est=Para_Post_pdf{i};
Summary = summarize(Est);
Summary = Summary.MarginalDistributions;
    Summary_new=Summary.Mean(1:end-2);
    Mean=Summary_new(:,1);
    Act_Inh_strength(i,setdiff(1:size(S,1),i))=Mean;
end
S_new=S;
for i=1:size(S,1)
    S_new(i,i)=0;
    S_new(i,setdiff(1:size(S,1),i))=S(i,1:end-1);
end
AM=Act_Inh_strength.*(S_new>0.95); 
csvwrite('Results\AdjacentMatrix_LPS.csv', AM); 