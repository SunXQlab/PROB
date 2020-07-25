%%%% Executable functions

%% Load Data
load ExampleData.mat  % Firstly load clinical cross-sectionla gene expression data;
Data = ExampleData; % The ExampleData is a dataset containing expression levels of 6 genes (the first 6 rows) and  clinical information (the last row).

%%  pseudotemporal progression inference
[Data_ordered,PPD,TimeSampled]=Progression_Inference(Data);

%% Bayesian Lasso for GRN parameter estimation 
[Para_Post_pdf,S,AM]=ODE_BayesianLasso(Data_ordered,TimeSampled);
csvwrite('AdjacentMatrix.csv', AM); % save Adjacent matrix for visualization (can be reformed for cytoscape visualization)

%% Visualization
Progression_Plot(X_stage);  % plot pseudotemporal progression trajectory

geneID=[1:6]; % specify the ID of genes you want to plot
TemporalGene_Plot(geneID,TimeSampled,Data_ordered); % plot pseudotemporal gene expression