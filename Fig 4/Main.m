%% Data load 

Expression=csvread('F:\Clinical Gene expression network Project\Reversion\Data\UC-SARC\GSE128192_Expression_processed.csv',1,1);
Staging=csvread('F:\Clinical Gene expression network Project\Reversion\Data\UC-SARC\Staging.csv',1,1);
EMT_reggene_ind=csvread('F:\Clinical Gene expression network Project\Reversion\Data\UC-SARC\EMT_reggene_ind.csv',1,1);
EMT_reggene_ind=EMT_reggene_ind(EMT_reggene_ind~=0);
EMT_reggene_ind=unique(EMT_reggene_ind,'stable');

input_data=[Expression;Staging'];
cd('F:\Clinical Gene expression network Project\Reversion\Codes\UC_SARC')
% [val,ind]=sort(Capture_time);
% input_data=input_data(:,ind);
save input_data.mat input_data
%% Progression Inferrence
tic
[Data_smooth,DPT,DPP]=Progression_Inferrence(input_data);
toc
% NetGene=Data_smooth;
[val,ind]=sort(DPT);
output_data=input_data(1:end-1,ind);
% output_data=Data_smooth(1:end-1,:);
save output_data.mat output_data
%%% 
%% selected TCGs
% Lcoeff=zeros(1,size(Data_smooth,1));
% Sdetrend=zeros(1,size(Data_smooth,1));
% for i=1:size(Data_smooth,1)
% p=polyfit(DPP,Data_smooth(i,:),1);
% Lcoeff(i)=p(1);
% Sdetrend(i)=std(detrend(Data_smooth(i,:)));
% end
% R=Lcoeff./Sdetrend; %Sdetrend;
% % R=1./Sdetrend; %Sdetrend;
% % R=Lcoeff./delta;
% [a,b]=sort(abs(R));
% ind_selected=b((end-99):end);
% Data_selected=output_data(ind_selected,:);
% % csvwrite('Results\ind_selected.csv', ind_selected);
% NetGene=Data_selected;
%% ODE_BayesianLasso network 
NetGene=output_data(EMT_reggene_ind,:);
DPT_sort=val';

 %% 
 tic
[Para_Post_pdf,S]=ODE_BayesianLasso(NetGene,DPT_sort);
toc
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
csvwrite('Results\AdjacentMatrix_UC_SARC.csv', AM); 

%% Differential network
NetGene1=NetGene(:,1:84);
NetGene2=NetGene(:,85:112);
DPT_sort1=DPT_sort(:,1:84);
DPT_sort2=DPT_sort(:,85:112);
save DPT_sort1.mat DPT_sort1 
save DPT_sort2.mat DPT_sort2
save NetGene1.mat NetGene1
save NetGene2.mat NetGene2
xlswrite('F:\E Disk\0-已完成的研究课题\PLoS CB - 2021\Manuscript\Codes\Results\Fig 4\NetGene1.xlsx',NetGene1)
xlswrite('F:\E Disk\0-已完成的研究课题\PLoS CB - 2021\Manuscript\Codes\Results\Fig 4\NetGene2.xlsx',NetGene2)
%%% 
[Para_Post_pdf1,S1]=ODE_BayesianLasso(NetGene1,DPT_sort1);
[Para_Post_pdf2,S2]=ODE_BayesianLasso(NetGene2,DPT_sort2);
%%% Adjacent matrix for visualization (cytoscape)
clear Mean_Para1  Mean_Para2
Act_Inh_strength_1=zeros(size(S1));Act_Inh_strength_2=zeros(size(S2));
for i=1:size(S1,1)
Est=Para_Post_pdf1{i};
Summary = summarize(Est);
Summary = Summary.MarginalDistributions;
    Summary_new=Summary.Mean(1:end-2);
    Mean=Summary_new(:,1);
    Act_Inh_strength_1(i,setdiff(1:size(S1,1),i))=Mean;
end

for i=1:size(S2,1)
Est=Para_Post_pdf2{i};
Summary = summarize(Est);
Summary = Summary.MarginalDistributions;
    Summary_new=Summary.Mean(1:end-2);
    Mean=Summary_new(:,1);
    Act_Inh_strength_2(i,setdiff(1:size(S2,1),i))=Mean;
end

S1_new=S1;
for i=1:size(S1,1)
    S1_new(i,i)=0;
    S1_new(i,setdiff(1:size(S1,1),i))=S1(i,1:end-1);
end

S2_new=S2;
for i=1:size(S2,1)
    S2_new(i,i)=0;
    S2_new(i,setdiff(1:size(S2,1),i))=S2(i,1:end-1);
end

AM1=Act_Inh_strength_1.*(S1_new>0.95); 
AM2=Act_Inh_strength_2.*(S2_new>0.95); 

csvwrite('Results\AdjacentMatrix_UC.csv', AM1); 
csvwrite('Results\AdjacentMatrix_SARC.csv', AM2); 

%% for trend visualization
Smoothed_data=Data_smooth(1:end-1,:);
NetGene_smoothed=Smoothed_data(EMT_reggene_ind,:);
save NetGene_smoothed.mat NetGene_smoothed
save DPP.mat DPP

