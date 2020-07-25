function Progression_Plot(X_stage)
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
for rep_root=1:size(data,1)   %randperm(size(data,1))
    a=sum(grade~=max(grade))+1;
    b=length(grade);
    rn=ceil(rand*(b-a)+a);
    drn=dpt_to_root(M,rn);
    [~,rr_123]=max(drn);  % select a node that has the maximal distance to a randomly selected node with maximal grade. 
    if grade(rr_123)~=max(grade)
    root=rr_123;
    break
    end
end
  
PPD=dpt_to_root(M,root)';

% compute diffusion components for plotting
[phi, lambda] = eig_decompose_normalized(T,5);

%plot
% label=zeros(1,length(PPD));
% label(1:5)=1;
% label(6:530+6)=2;
% label(530+7:length(PPD))=3;

% label=grade;
% 
figure,
% scatter3(phi(:,2),phi(:,3),phi(:,4),50,PPD,'fill');
scatter(-phi(:,2),phi(:,3),50,PPD,'fill');
% scatter(-phi(:,2),phi(:,3),50,label,'fill');
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('DC1')
ylabel('DC2')
set(gca,'FontSize',20)
colorbar
plot_average_path(PPD,1:length(PPD),[-phi(:,2),phi(:,3)])
% % 
figure,
plot(1:length(PPD),sort(PPD),'-','color',[0.5 0.5 0.5],'LineWidth',2)
hold on
scatter(1:length(PPD),sort(PPD),50,sort(PPD),'fill')
set(gca,'FontSize',15)
xlabel('Order in pseudo-progression','FontSize',20);
ylabel('Pseudo-progression score','FontSize',20)

