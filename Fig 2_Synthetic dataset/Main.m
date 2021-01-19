%%%% Fig 2-associated code

ODESimulator;
[Data_smooth,DPT,DPP]=Progression_Inferrence(X_stage);
[Para_Post_pdf,S]=ODE_BayesianLasso(Data_smooth,DPP);
   for i=1:size(C,1)
        S_new(i,i)=0;
        S_new(i,setdiff(1:size(X,1),i))=S(i,1:end-1);
   end
    
   figure,
   ROC_data4=plot_roc(S_new, C~=0 ,1,1);  % ROC calculation for PPBayLasso
