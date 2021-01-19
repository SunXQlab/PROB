
%%% Main function
cd('H:\E Disk\Clinical Gene expression network Project\Manuscript\Codes\Fig 4-5_Realistic applications\Breast cancer TCGA_Fig 4\')
%% BRCA dataset
Gene=xlsread('Data\Gene_TCGA_BRCA_processed.csv');
Stage=xlsread('Data\Phenotype_TCGA_BRCA_processed.csv');
size(Gene)
size(Stage)
X_stage=[Gene;Stage';1:size(Stage,1)];
size(X_stage)
%% Pseduo_Progression
[Data_smooth,DPT,DPP]=Progression_Inference(X_stage);
 
%% selected genes ordered by DPT
Gene_Name=textread('Breast cancer\SelectedGeneName_GEO_BreastCancer.txt','%s');  %,'whitespace',''
Gene_Index=xlsread('Breast cancer\Gene_selected_Index.csv');  %,'whitespace',''
Gene_Index=Gene_Index(:,2);
Data_selected=Data_smooth(Gene_Index,:);

csvwrite('Breast cancer\GeneExpression_Selected_Smoothed.csv', Data_selected); 
csvwrite('Breast cancer\GeneExpression_Selected_NonSmoothed.csv', Gene(Gene_Index,:)); 
%% Plot time-couse 
[value,ind]=sort(DPT);
figure,
heatmap((zscore((Data_selected(:,1:50:1000))'))');
colormap cool
snapnow

% GID = R_ind(1:10);
% 
% for i=1:length(GID)  
% figure,
% smoothL=100;
% g_smooth=ksmooth(data(indT,GID(i)),smoothL);
% % [g_smooth,window]=smoothdata(data(:,g),'gaussian');
% scatter(DPP,data(indT(ceil(smoothL/2):ceil(length(DPT)-smoothL/2)-1),GID(i)),20,[0.5 0.5 0.5],'fill')
%  hold on
% 
% Time_smooth=valT(ceil(smoothL/2):ceil(length(DPT)-smoothL/2)-1);  %CEIL(X) rounds the elements of X to the nearest integers towards infinity.
% DPP=(ceil(smoothL/2):ceil(length(DPT)-smoothL/2)-1)/max(ceil(smoothL/2):ceil(length(DPT)-smoothL/2)-1);
% scatter(DPP,g_smooth,50,Time_smooth,'fill')
% hold on
% % colormap bone;brighten(0.6)
% plot(DPP,g_smooth,'-','color',[0.5 0.5 0.5],'LineWidth',1)  %  
% set(gca,'FontSize',15)
% xlabel('Order in pseudo-progression','FontSize',20);
% ylabel('Gene expression','FontSize',20)
% % title(G_Name(i))
% end

figure, 
plot(DPP,Data_selected(71,:))

%% MI Network
N=size(Data_selected,1);
MI=zeros(N,N);
Gene_selected=Gene(Gene_Index,:);
for i=1:N
    for j=1:N
        MI(i,j)=VectorMI(Gene_selected(i,:)',Gene_selected(j,:)');
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

csvwrite('Breast cancer\InteractionScore.csv', S_new);  
save ('Breast cancer\Para_Post_pdf.mat','Para_Post_pdf')

for i=1:72
Act_Inh_strength(i,:)=(Para_Post_pdf{i}.Mean)';
end

%% Adjacent matrix for visualization (cytoscape)
AM=Act_Inh_strength.*(S_new>0.95); 
csvwrite('Breast cancer\AdjacentMatrix.csv', AM); 


%%%

corrcoef(Gene_selected(7,:)',Gene_selected(71,:)')


