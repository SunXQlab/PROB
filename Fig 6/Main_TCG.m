
%%% Main function
cd('H:\E Disk\Clinical Gene expression network Project\Manuscript\Codes\Fig 4-5_Realistic applications\Breast cancer GEO_Fig 5\')
%% BRCA dataset
Gene=xlsread('Data\Gene_GSE7390.csv');
Grade=xlsread('Data\Grade_GSE7390.csv');
% Gene=Gene(:,2:end);
size(Gene)
size(Grade)
X_Grade=[Gene;Grade';1:size(Grade,1)];
size(X_Grade)
%% Pseduo_Progression
[Data_smooth,DPT,DPP]=Progression_Inference(X_Grade);
 
%% selected TCGs
Lcoeff=zeros(1,size(Data_smooth,1));
Sdetrend=zeros(1,size(Data_smooth,1));
for i=1:size(Data_smooth,1)
p=polyfit(DPP,Data_smooth(i,:),1);
Lcoeff(i)=p(1);
Sdetrend(i)=std(detrend(Data_smooth(i,:)));
end
R=Lcoeff./Sdetrend; %Sdetrend;
% R=1./Sdetrend; %Sdetrend;
% R=Lcoeff./delta;
[a,b]=sort(abs(R));
ind_selected=b((end-99):end);
Data_selected=Data_smooth(ind_selected,:);

figure,plot(1:length(value),Gene(ind_selected(end),ind))
csvwrite('Results\ind_selected.csv', ind_selected);

csvwrite('Results\Data_smooth_selected.csv', Data_selected);


%% MI Network
N=size(Data_selected,1);
MI=zeros(N,N);
for i=1:N
    for j=1:N
        MI(i,j)=VectorMI(Data_selected(i,:)',Data_selected(j,:)');
    end
end
MI(logical(eye(size(MI,1))))=0;
Threshold=sort(reshape(MI,1,N*N),'descend');
MI=(MI>=Threshold(floor(N*N*0.05)));
% MI(logical(eye(size(MI,1))))=1;

%% ODE_BayesianLasso
% cd('F:\Clinical Gene expression network Project\Code\TCGA Data\')
% PriorPPI=xlsread('Breast cancer\MI Network.csv');
tic
[Para_Post_pdf,S]=ODE_BayesianLasso_PriorPPI(Data_selected,DPP,MI);
toc

%% Adjusted interaction matrix
S_new=S;
% for i=1:size(S,1)
%     S_new(i,i)=0;
%     S_new(i,setdiff(1:size(S,1),i))=S(i,1:end-1);
% end
% 
csvwrite('Results\InteractionScore.csv', S_new);  
save ('Results\Para_Post_pdf.mat','Para_Post_pdf')

%%  SVD for measuring node importance
clear Mean_Para Var_Para
for i=1:size(S_new,1)
Mean_Para(i,:)=table2array(Para_Post_pdf{i}(:,1));  % mean 
Var_Para(i,:)=table2array(Para_Post_pdf{i}(:,2));  % variation 
end

[X,Y]=eig(Mean_Para*Mean_Para')
diag(Y)
X(:,end)
[a,b]=sort(X(:,end))

ind_selected(b(end))

csvwrite('Results\Node_Importance.csv', [a,ind_selected(b)']); 
%% Adjacent matrix for visualization (cytoscape)
clear Act_Inh_strength
for i=1:size(S_new,1)
Act_Inh_strength(i,:)=(Para_Post_pdf{i}.Mean)';
end

AM=Act_Inh_strength.*(S_new>0.95); 
csvwrite('Results\AdjacentMatrix.csv', AM); 

%% Plot time-couse 
DDD=(Data_smooth(ind_selected,:)-min(Data_smooth(ind_selected,:),[],2))./(max(Data_smooth(ind_selected,:),[],2)-min(Data_smooth(ind_selected,:),[],2));
figure, plot(DPP,DDD)

DDD=(Gene(ind_selected,:)-min(Gene(ind_selected,:),[],2))./(max(Gene(ind_selected,:),[],2)-min(Gene(ind_selected,:),[],2));
figure, plot(1:196,DDD)

[gg,pp]=sort(Grade);
figure, plot(Gene(ind_selected(b(end)),pp))

[gg,pp]=sort(Grade);
figure, plot(Gene(ind_selected,pp))

%% 
