function [Para_Post_pdf,S,AM]=ODE_BayesianLasso(Data_ordered,TimeSampled)
% Input:
% Data_smooth -- (pseudo)time-course gene expression data
% TimeSampled -- Time points associate with (pseudo)time-course gene expression data
% Output:
% Para_Post_pdf -- posterior distribution over the coefficients in the GRN model
% S --  a matrix saving the presence probability. 
% AM -- Adjacent matrix of the inferred network i.e., (a_ij) for the regulatory coefficient from gene j to gene i. 
%% data reformation
Time=TimeSampled;
x=Data_ordered;  % original data 
Time=(Time-min(Time))/(max(Time)-min(Time)); % standardize
y=diff((x)')./diff(Time)';  % derivative term approximation
y=y';
x=x(:,1:end-1);
%% Calculate Post_pdf for each gene's parameters
Para_Post_pdf=cell(size(x,1));
for ind_gene=1:size(x,1)
    
y_output=y(ind_gene,:);
x_input=[x',ones(size(x,2),1)].*x(ind_gene,:)';
x_input(:,ind_gene)=[];    
%% Lasso
% [A,FitInfo] = lasso_Sun(x_input,y_output,'Alpha',1,'CV',10,'Standardize',1,'Intercept_Flag',0);  % modified from matlab function by adding Intercept_Flag argument associated with lassoFit_Sun (redefined function)
% idxLambdaMinMSE = FitInfo.IndexMinMSE;  %Predict exam scores for the test data. 
% coef = A(:,idxLambdaMinMSE);
% coef0 = FitInfo.Intercept(idxLambdaMinMSE);
%% Bayesian Lasso
% predictornames={'Deg','B','C'};
p=size(x,1);
PriorMdl = bayeslm(p,'ModelType','lasso','Intercept',false);  %,'Lambda',L  %% ,'VarNames',predictornames; %% ,'Mu',[coef],'V',eye(p)
% table(PriorMdl.Lambda,'RowNames',PriorMdl.VarNames)
ismissing = any(isnan(x),2);
n = sum(~ismissing); % Effective sample size

% Preallocate
% BayesLassoCoefficients = zeros(p,1);
% BayesLassoCI95 = zeros(p,2,1);
% fmseBayesLasso = zeros(1,1);
% BLCPlot = zeros(p,1);

% Estimate and forecast
rng(3); % For reproducibility
Y_train=y_output;
X_train=x_input;
   
[EstMdl,Summary] = estimate(PriorMdl,X_train,Y_train,'Display',false);  %,'BetaStart',[coef],'Sampler',"hmc",'NumDraws',1e6,'Burnin',1e3

% figure, plot(EstMdl)  %%%% Plot the posterior distribution of estimated  parameters
Para_Post_pdf{ind_gene}=EstMdl;   % save the parameter estimation for targeting each gene into a cell structure
end

%%  Calculate Confidence Interval
    alpha=0.01:0.01:1;
    CI_alpha=zeros(size(x,1),length(alpha),size(x,1)+1,2);
for i=1:size(x,1)
    Para_Post_pdf_each=Para_Post_pdf{i};
%     Para=Para_Post_pdf_each.BetaDraws';
Summary = summarize(Para_Post_pdf_each);
Summary = Summary.MarginalDistributions;
    alpha=0.01:0.01:1;
    ind_alpha=1;
 for k_alpha=alpha
%     MUCI = CI_Sun(Para,k_alpha);  % 
    MUCI = CI_Sun_Normal(Summary,k_alpha);  % 
    CI_alpha(i,ind_alpha,:,:)=MUCI;
    ind_alpha=ind_alpha+1;
 end
    size(CI_alpha);
end

S=ones(size(x,1),size(x,1));  % the presence probability of each predicted interaction; Intercept term counted
for i=1:size(x,1) 
    for k=length(alpha):-1:1
        CI_genei_alpha=reshape(CI_alpha(i,k,:,:),size(x,1)+1,2);
        for j=1:size(x,1)
             if CI_genei_alpha(j,1)<=0 & CI_genei_alpha(j,2)>=0   % for coeficient matric               
                 S(i,j)=min(S(i,j),1-alpha(k));  % S=1-max{alpha:CI_alpha containing 0}
%                  break;
             end        
       end

   end
end

%% Reform probability matrix for interation coefficients
   for i=1:size(Data_ordered,1)
        S_new(i,i)=0;
        S_new(i,setdiff(1:size(Data_ordered,1),i))=S(i,1:end-1);
   end
   S=S_new;

   
 %% Output for network visualization
clear Act_Inh Act_Inh_strength
for i=1:size(S,1)
Summary = summarize(Para_Post_pdf{i});  %(Para_Post_pdf{i}.Mean)';
Summary = Summary.MarginalDistributions;
Act_Inh(i,:) = Summary.Mean(1:end-2);
Act_Inh_strength(i,i)=0;
Act_Inh_strength(i,setdiff(1:size(S,1),i))=Act_Inh(i,:);
end
% Adjacent matrix
AM=Act_Inh_strength.*(S>0.95);