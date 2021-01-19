
function [Data_ordered,PPD,TimeSampled]=Progression_Inference(X_stage)
%%
% X_stage is clinical cross-sectional dataset containing gene expression profiles and grade information of patients.
% Data_smooth is ordered and smoothed gene expression data along with pseudotemporal progression trajectory. 
% PPD is pseudotemporal progression distance for each patient with respect to the rooting sample.
% TimeSampled is standardized time-points sampled for Data_smooth.
%%
R=size(X_stage);
data=X_stage(1:R(1)-1,:)';
% data=cosine_dist_normalisation(data);
grade=X_stage(R(1),:);

%%  clinical information-weighted locally-scaled kernel
% weight coefficients
W=zeros(length(grade),length(grade));
for i=1:length(grade)
    for j=1:length(grade)
    W(i,j)=(1+abs(grade(i)-grade(j)));
    end
end

% make transition matrix
nsig=10; [T,phi0] = T_loc(data,nsig,W);

%make the accumulated transition matrix and find tip cells
M = dpt_input(T, phi0);

%% calculate root
% for rep_root=1:size(data,1)   %randperm(size(data,1))
%     a=sum(grade~=max(grade))+1;
%     b=length(grade);
%     rn=ceil(rand*(b-a)+a);
%     drn=dpt_to_root(M,rn);
% 
%     [~,rr_123]=max(drn);  % select a node that has the maximal distance to a random selected node in grade 4. 
%     if grade(rr_123)==min(grade) %~=4
%     root=rr_123;
% %     rep_root
%     break
%     end
% end
 %% another approach to select root
    Ind_max=find(grade==max(grade));
    b=length(max(grade));
    rn=ceil(rand*b);
    x_ref=Ind_max(rn);
    drn=dpt_to_root(M,x_ref);
    
    [~,AA_sort]=sort(drn,'descend');
    AA=find(grade(AA_sort)==min(grade));   % limit the selection of root to Grade_min
    root=AA(1);
 
  
PPD=dpt_to_root(M,root)';


[valT,indT]=sort(PPD);
smoothL=10^(floor(log10(R(2)))-1);
g=1; 
g_smooth=ksmooth(data(indT,g),smoothL);
Time_smooth=valT(ceil(smoothL/2):ceil(length(PPD)-smoothL/2)-1);  %CEIL(X) rounds the elements of X to the nearest integers towards infinity.


Data_ordered=zeros(size(g_smooth,1),size(data,2));
for i=1:size(data,2)
Data_ordered(:,i)=ksmooth(data(indT,i),smoothL);
end
Data_ordered=Data_ordered';

TimeSampled=[smoothL/2+1:size(X_stage,2)-smoothL/2]/size(X_stage,2);  % Diffusion pseudo-progression  % size(X_stage,2)=100
% TimeSampled=0:1e-2:1;
% TimeSampled=1:size(Data_smooth,2);
Data_ordered=reshape(Data_ordered,size(data,2),length(TimeSampled));

