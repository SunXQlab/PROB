%%%% Fig 2-associated code

scExpression=csvread('F:\Clinical Gene expression network Project\Reversion\Data\LPS_WT_scExpression.csv',1,1);
Capture_time=csvread('F:\Clinical Gene expression network Project\Reversion\Data\LPS_Capture_time.csv',1,1);

Regulators_Targets_ind=csvread('F:\Clinical Gene expression network Project\Reversion\Data\Regulators_Targets_ind.csv',1,1);

input_data=[scExpression;Capture_time'];
cd('F:\Clinical Gene expression network Project\Reversion\Codes\Benchmarking_LPS_scRNAseq_data')
save scExpression_Time.mat input_data
[val,ind]=sort(Capture_time);
input_data=input_data(:,ind);
[Data_smooth,DPT,DPP]=Progression_Inferrence(input_data);
Regulators_Targets=Data_smooth(Regulators_Targets_ind,:);
[Para_Post_pdf,S]=ODE_BayesianLasso(Regulators_Targets,DPP);
   save DPT.mat DPT
   scExp_Net_Ctime=[input_data(Regulators_Targets_ind,:);DPT'];
   save scExp_Net_Ctime.mat scExp_Net_Ctime
   
   for i=1:size(Regulators_Targets,1)
        S_PROB_RT(i,i)=0;
        S_PROB_RT(i,setdiff(1:size(Regulators_Targets,1),i))=S(i,1:end-1);
   end
  save S_PROB_RT.mat  S_PROB_RT
  
clear Act_Inh Act_Inh_strength
for i=1:size(S,1)
Summary = summarize(Para_Post_pdf{i});  %(Para_Post_pdf{i}.Mean)';
Summary = Summary.MarginalDistributions;
Act_Inh(i,:) = Summary.Mean(1:end-2);
Act_Inh_strength(i,i)=0;
Act_Inh_strength(i,setdiff(1:size(S,1),i))=Act_Inh(i,:);
end
% Adjacent matrix
AM=Act_Inh_strength.*(S_PROB_RT>0.95);
save AM.mat AM
save Act_Inh_strength.mat Act_Inh_strength

