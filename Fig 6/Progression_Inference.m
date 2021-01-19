
function [Data_smooth,DPT,DPP]=Progression_Inference(X_stage)
% X_stage=xlsread('F:\Clinical Gene expression network Project\Code\Simulated ODE data\GeneExpression_Stage_Simulated.csv');
% load X_stage.mat
% clear data
% cd ('F:\Clinical Gene expression network Project\Code\Simulated ODE data - 3GeneNet\Noise_Effect')
R=size(X_stage);
data=X_stage(1:R(1)-2,:)';
% data=cosine_dist_normalisation(data);
grade=X_stage(R(1)-1,:);
% grade=X_stage(R(1),:);  % real order
W=zeros(length(grade),length(grade));
for i=1:length(grade)
    for j=1:length(grade)
    W(i,j)=(1+abs(grade(i)-grade(j)));
%     W(i,j)=1;
%     W(i,j)=0.1*W(i,j);
    end
end

% root=find(rand_index==1); %specify the index of the root cell

% make transition matrix
%k=50; nsig=10; [T,phi0] = diffusionmap.T_nn(data,data,k,nsig);

nsig=3; [T,phi0] = T_loc(data,nsig,W);
% sigma=0.0005; [T,phi0] = T_classic(data,sigma);
%make the accumulated transition matrix and find tip cells
M = dpt_input(T, phi0);


% calculate root
% clear root
for rep_root=1:size(data,1)   %randperm(size(data,1))
    a=sum(grade~=max(grade))+1;
    b=length(grade);
    rn=ceil(rand*(b-a)+a);
    drn=dpt_to_root(M,rn);
    [~,rr_123]=max(drn);  % select a node that has the maximal distance to a random selected node in grade 4. 
    if grade(rr_123)~=max(grade)
    root=rr_123;
    break
    end
end
  
DPT=dpt_to_root(M,root)';

% compute diffusion components for plotting
[phi, lambda] = eig_decompose_normalized(T,5);


%%%%%%%%%%%% plot expression of genes along the DPT
[valT,indT]=sort(DPT);
% figure,
smoothL=10^(floor(log10(R(2)))-1);
g=1; %19387; %18206; %3318;
g_smooth=ksmooth(data(indT,g),smoothL);
% [g_smooth,window]=smoothdata(data(:,g),'gaussian');

% scatter(1:length(DPT),data(:,g),20,'b','fill')
%  hold on

Time_smooth=valT(ceil(smoothL/2):ceil(length(DPT)-smoothL/2)-1);  %CEIL(X) rounds the elements of X to the nearest integers towards infinity.



%%%%%%%%%%%%  Save the DPT and the smoothed data
% csvwrite('F:\Clinical Gene expression network Project\Code\Simulated ODE data\DPT_smooth.txt', Time_smooth);  
% save Time_smooth.mat

Data_smooth=zeros(size(g_smooth,1),size(data,2));
for i=1:size(data,2)
Data_smooth(:,i)=ksmooth(data(indT,i),smoothL);
end
Data_smooth=Data_smooth';
% DPP=[smoothL/2+1:size(X_stage,2)-smoothL/2]/size(X_stage,2);  % Diffusion pseudo-progression  % size(X_stage,2)=100
% DPP_new=min(DPP):1e-2:1;
% Data_smooth_Interpolate=pchip(DPP,Data_smooth',DPP_new);  %% 三次Hermite多项式插值，可用pchip函数替代
% clear Data_smooth_Inter_Extrap
% for i=1:size(Data_smooth',1)
% Data_smooth_extrap=interp1(DPP_new(3:end),Data_smooth_Interpolate(i,3:end)',[0:1e-2:(min(DPP)+1*1e-2)],'linear','extrap');  %% ,'extrap'
% Data_smooth_Inter_Extrap(i,:)=[Data_smooth_extrap,Data_smooth_Interpolate(i,3:end)];
% end
% 
% Data_smooth=Data_smooth_Inter_Extrap;


% csvwrite('F:\Clinical Gene expression network Project\Code\Simulated ODE dataGeneExpression_DPT_smooth.txt', Data_smooth);   
% save Data_smooth.mat
% SS=writetable( Data_smooth,'F:\Clinical Gene expression network Project\Results\Data_smooth.csv');

%%%%%%%%%%%  Plot time course curves

DPP=[smoothL/2+1:size(X_stage,2)-smoothL/2]/size(X_stage,2);  % Diffusion pseudo-progression  % size(X_stage,2)=100
% DPP=0:1e-2:1;
% DPP=1:size(Data_smooth,2);
Data_smooth=reshape(Data_smooth,size(data,2),length(DPP));

