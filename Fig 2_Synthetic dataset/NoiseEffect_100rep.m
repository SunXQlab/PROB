close all
clear all
clc

tic 
global NoiseLevel
 
grid_noiselevel=0:1:30;
ind=1;
for NoiseLevel=grid_noiselevel  %  x%

    NoiseLevel
    for rep=1:100   %% repeat random simulation 100 times
        rep
    ODESimulator;
    [Data_smooth,DPT,DPP]=Diffusion_Pseudo_Progression(X_stage);
    [Para_Post_pdf,S]=ODE_BayesianLasso(Data_smooth,DPP);
    
%% ROC calculation
    S_new=S;
    for i=1:size(C,1)
        S_new(i,i)=0;
        S_new(i,setdiff(1:size(X,1),i))=S(i,1:end-1);
    end
    C(logical(eye(size(C,1))))=0;
    ROC_data=plot_roc(S_new, C~=0 ,0,0);  % [C,D'] correspond to Intercept term and C 
    
%% Evaluation Index
    
       
    TPR_noise(ind,rep)=ROC_data.param.Sensi;  % Sensi is TPR
    FPR_noise(ind,rep)=ROC_data.param.FPR;
    AUC_noise(ind,rep)=ROC_data.param.AROC;
    Accuracy_noise(ind,rep)=ROC_data.param.Accuracy;
    PPV_noise(ind,rep)=ROC_data.param.PPV;
    MCC_noise(ind,rep)=ROC_data.param.MCC;
    S_distance(ind,rep)=sqrt(sum((X_stage(8,:)/100-(DPT'/max(DPT'))).^2)/length(DPT'));
    coeff_Spearman(ind,rep) = corr(X_stage(8,:)', DPT , 'type' , 'Spearman');
    coeff_Kendall(ind,rep) = corr(X_stage(8,:)', DPT , 'type' , 'Kendall');
    end
    ind=ind+1;
end
toc

figure,
errorbar(grid_noiselevel,mean(AUC_noise,2),std(AUC_noise,0,2),'o-.','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.3 0.8 0.5]);title('AUC'); axis([-1 31 0 1]); 

figure,
errorbar(grid_noiselevel,mean(TPR_noise,2),std(TPR_noise,0,2),'o-.','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.3 0.8 0.5]); title('TPR'); axis([-1 31 0 1]);
figure,
errorbar(grid_noiselevel,mean(FPR_noise,2),std(FPR_noise,0,2),'o-.','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.3 0.8 0.5]); title('FPR'); axis([-1 31 0 1]);
figure,
errorbar(grid_noiselevel,mean(Accuracy_noise,2),std(Accuracy_noise,0,2),'o-.','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.3 0.8 0.5]); title('Accuracy'); axis([-1 31 0 1]);
figure,
errorbar(grid_noiselevel,mean(PPV_noise,2),std(PPV_noise,0,2),'o-.','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.3 0.8 0.5]); title('PPV'); axis([-1 31 0 1]); % axis([-1 31  0 max(mean(PPV_noise,2)+std(PPV_noise,0,2))])
figure,
errorbar(grid_noiselevel,mean(MCC_noise,2),std(MCC_noise,0,2),'o-.','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.3 0.8 0.5]); title('MCC'); axis([-1 31 0 1]); % axis([-1 31  0 max(mean(MCC_noise,2)+std(MCC_noise,0,2))])

figure,
errorbar(grid_noiselevel,mean(S_distance,2),std(S_distance,0,2),'o-.','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.1 0.3 1]);title('RMSE of Progression Inference'); axis([-1 31 0 0.3]); % axis([-1 31  0 max(mean(S_distance,2)+std(S_distance,0,2))])

figure,
errorbar(grid_noiselevel,mean(coeff_Spearman,2),std(coeff_Spearman,0,2),'o-.','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.1 0.3 1]);title('Spearman correlation'); axis([-1 31 0.5 1]); % axis([-1 31  0 max(mean(S_distance,2)+std(S_distance,0,2))])

figure,
errorbar(grid_noiselevel,mean(coeff_Kendall,2),std(coeff_Kendall,0,2),'o','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.1 0.3 1]);title('Kendall correlation of Progression Inference'); axis([-1 31 0 1]); % axis([-1 31  0 max(mean(S_distance,2)+std(S_distance,0,2))])

%% setting for each figure
set(gca,'FontSize',15)
xlabel('Noise level (%)','FontSize',20);

% figure,
% boxplot(S_distance','plotstyle','compact','color',[0.5 0.5 1]);axis([-1 51  0 max(max(S_distance))])
% set(gca,'XTickLabel',sprintfc('%g',0:5:50))



% 
% 
% figure,
% plot(0:5:50,AUC_noise); title('AUC'); axis([0 50 0 1])
% 
% figure,
% plot(0:5:50,TPR_noise); title('TPR'); axis([0 50 0 1])
% 
% figure,
% plot(0:5:50,FPR_noise); title('FPR'); axis([0 50 0 1])
% 
% figure,
% plot(0:5:50,Accuracy_noise); title('Accuracy'); axis([0 50 0 1])
% 
% figure,
% plot(0:5:50,PPV_noise); title('PPV'); axis([0 50 0 1])
% 
% figure,
% plot(0:5:50,MCC_noise); title('MCC'); axis([0 50 0 1])
% 
%     
%     