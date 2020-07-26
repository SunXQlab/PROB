# PROB
pseudotemporal progression-based Bayesian method for inferring GRNs from cross-sectional clinical transcriptomic data

PROB contains two main functions, ‘Progression_Inferrence’ (for inferring pseudotemporal progression) and ‘ODE_BayesianLasso’ (for inferring ODE network using Bayessian Lasso). The visualization module of PROB enables the users to plot pseudotemporal progression trajectory and to plot pseudotemporal expression of the selected genes. In addition, the prior network information can be incorporated into function ‘ODE_BayesianLasso’ when PROB is applied to large scale network reconstruction. 

The input and output of each function are described below. 

## 1. Pseudotemporal progression inference

**[Data_ordered,PPD, TimeSampled]=Progression_Inferrence(Data)**

This function (.m) is designed to infer pseudotemporal progression from the clinical cross-sectional gene expression data. 

***Input***

Data: a n×m matrix containing gene expression profiles (the first n-1 rows) and grade information of patients (the last row). m is sample size, i.e., the number of patients. 

***Output***

Data_ordered: a matrix containing the ordered (and smoothed) gene expression data.

PPD: a vector containing pseudotemporal progression distance for each patient.

TimeSampled: a vector containing standardized time-points sampled for Data_ordered. 

## 2. GRN inference

**[Para_Post_pdf,S,AM]=ODE_BayesianLasso(Data_ordered, TimeSampled)**

This function (.m) is designed to reconstruct causal GRN based on the results of the above pseudotemporal progression inference. 

***Input***

Data_ordered: matrix for the ordered (and smoothed) gene expression data (i.e., a subset from the output of the first function).

TimeSampled: vector for the standardized time-points associated with Data_ordered (i.e., the output of the first function).  

***Output***

Para_Post_pdf: cell format saving the posterior distribution over the regulatory coefficients of each gene in the GRN model. 

S: a matrix saving the presence probability. 

AM: Adjacent matrix of the inferred GRN. (aij) for the regulatory strength from gene j to gene i. 

## 3. Visualization

* **Progression_Plot(Data)**

This function (.m) is designed to plot pseudotemporal progression trajectory. 

***Input***

Data: a n×m matrix containing single cell gene expression profiles (the first n-1 rows) and grade information of patients (the last row). m is sample size, i.e., the number of patients.

***Output***

Two figures of pseudotemporal progression trajectory and pseudotemporal progression score along the ordering of patients. 

* **TemporalGene_Plot(geneID,TimeSampled,Data_ordered)**

This function (.m) is designed to plot pseudotemporal expression of the selected genes.

***Input***

geneID: the ID of genes selected for visualization.

Data_ordered: matrix for the ordered and smoothed gene expression data (i.e., the output of the first function).

TimeSampled: vector for the standardized time-points associated with Data_ordered (i.e., the output of the first function).  

***Output***

A figure of the expression dynamics of the selected genes over the inferred pseudotemporal progression.

* **Cytoscape_Reformat(AM,NodeID)**

This function (.R) is designed to reformat the AM to a suitable format for input as cytoscape software. 

***Input***

AM: The adjacent matrix of the inferred GRN resulted from the above function ‘ODE_BayesianLasso’.

NodeID: The ID or the symbol of the genes in the AM.  

***Output***

A matrix containing 3 columns: source nodes, interaction coefficients and target nodes. 

## 4. Incorporating prior network information

**[Para_Post_pdf,S,AM]=ODE_BayesianLasso_PriorNet(Data_ordered,TimeSampled,PriorNet)**

This function (.m) is designed to incorporate prior network information into function ‘ODE_BayesianLasso’, when PROB is applied to large scale network reconstruction. 

***Input***

Data_ordered: matrix for the ordered and ordereded gene expression data (i.e., the output of the first function).

TimeSampled: vector for the standardized time-points associated with Data_ordered (i.e., the output of the first function).  

PriorNet: a matrix for the prior network information. (Pij) represents the regulatory strength from gene j to gene i.  

***Output***

Para_Post_pdf: cell format saving the posterior distribution over the regulatory coefficients of each gene in the GRN model. 

S: a matrix saving the presence probability. 

AM: Adjacent matrix of the inferred GRN. (aij) for the regulatory strength from gene j to gene i. 

## Example demonstration

The users need to first download all the files in this repository and save them into the working directory.

 *Load data*
    
    load ExampleData.mat  % Firstly load example data;
    Data = ExampleData; % The ExampleData is a dataset containing expression levels of 6 genes (the first 6 rows) and  clinical information (the last row).

 *pseudotemporal progression inference*
    
    [Data_ordered,PPD,TimeSampled]=Progression_Inference(Data);

 *Bayesian Lasso for GRN parameter estimation*
    
    [Para_Post_pdf,S,AM]=ODE_BayesianLasso(Data_ordered,TimeSampled);
    csvwrite('AdjacentMatrix.csv', AM); % save Adjacent matrix for visualization (can be reformed for cytoscape visualization)

 *Visualization*
    
    Progression_Plot(X_stage);  % plot pseudotemporal progression trajectory
    geneID=[1:6]; % specify the ID of genes you want to plot
    TemporalGene_Plot(geneID,TimeSampled,Data_ordered); % plot pseudotemporal gene expression

## 6. Other information

Please cite oue paper if you used codes here. 

Any questions please contact sunxq6@mail.sysu.edu.cn
