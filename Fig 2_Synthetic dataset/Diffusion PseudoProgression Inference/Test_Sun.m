
% Test

% data=load ('F:\Clinical Gene expression network Project\DPT_inMatlab\dpt\examples\ESC_qpcr_Goettgens.mat');
% d=struct2cell(data);
% data=d{3};

data1=xlsread('F:\Clinical Gene expression network Project\Data\Gene_Normal.xls');
data2=csvread('F:\Files from office\胶质瘤诱导分化组学数据分析\Code-cluster_based\Survival analysis\Data\Gene_Tumor.csv',1,1);
data=[data1,data2]';
LOD=-14; k=20; data=lam(data,LOD,k); %preprocessing
data=cosine_dist_normalisation(data);

grade_LGG=xlsread('F:\Clinical Gene expression network Project\Data\grade_LGG.xls');
grade_LGG= grade_LGG(:,2);
grade=[ones(1,size(data1,2)),grade_LGG',ones(1,size(data2,2)-530)*4];
W=zeros(length(grade),length(grade));
for i=1:length(grade)
    for j=1:length(grade)
    W(i,j)=(1+abs(grade(i)-grade(j)));
    end
end

root=1:5; %specify the index of the root cell

% make transition matrix
%k=50; nsig=10; [T,phi0] = diffusionmap.T_nn(data,data,k,nsig);
nsig=20; [T,phi0] = T_loc(data,nsig,W);
% sigma=0.05; [T,phi0] = T_classic(data,sigma);
%make the accumulated transition matrix and find tip cells
M = dpt_input(T, phi0);
% do pseudotimeordering and branch separation
% branching=1;
% [Branch,DPT]=dpt_analyse(M,branching,tips);

DPT=dpt_to_root(M,root)';

% compute diffusion components for plotting
[phi, lambda] = eig_decompose_normalized(T,5);

%plot
label=zeros(1,length(DPT));
label(1:5)=1;
label(6:530+6)=2;
label(530+7:length(DPT))=3;

label=grade;

figure,
% scatter3(phi(:,2),phi(:,3),phi(:,4),30,DPT,'fill');
% scatter(-phi(:,2),phi(:,3),50,DPT,'fill');
scatter(-phi(:,2),phi(:,3),50,label,'fill');
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('DC1')
ylabel('DC2')
set(gca,'FontSize',20)
colorbar
plot_average_path(DPT,1:length(DPT),[-phi(:,2),phi(:,3)])

figure,
plot(1:length(DPT),sort(DPT),'-','color',[0.5 0.5 0.5],'LineWidth',2)
hold on
scatter(1:length(DPT),sort(DPT),50,sort(DPT),'fill')
set(gca,'FontSize',15)
xlabel('Order in pseudo-progression','FontSize',20);
ylabel('Pseudo-progression score','FontSize',20)

%%%%%%%%%%%% plot expression of genes along the DPT
figure,
smoothL=200;
g=18206; %19387; %18206; %3318;
g_smooth=ksmooth(data(:,g),smoothL);
% [g_smooth,window]=smoothdata(data(:,g),'gaussian');

% scatter(1:length(DPT),data(:,g),20,'b','fill')
%  hold on
[valT,indT]=sort(DPT);
Time_smooth=valT(ceil(smoothL/2):ceil(length(DPT)-smoothL/2)-1);  %CEIL(X) rounds the elements of X to the nearest integers towards infinity.
scatter(Time_smooth,g_smooth/g_smooth(1),50,Time_smooth,'fill')
hold on
% colormap bone;brighten(0.6)
plot(Time_smooth,g_smooth/g_smooth(1),'-','color',[0.5 0.5 0.5],'LineWidth',1)  %  
set(gca,'FontSize',15)
xlabel('Order in pseudo-progression','FontSize',20);
ylabel('Gene expression','FontSize',20)
% ylim([min(data(:,g)) max(data(:,g))])


%%%%%%%%%%%%  Save the DPT and the smoothed data
xlswrite('F:\Clinical Gene expression network Project\Results\Time_smooth.xls', Time_smooth);  

tic
Data_smooth=zeros(size(g_smooth,1),size(data,2));
for i=1:size(data,2)
Data_smooth(:,i)=ksmooth(data(:,i),smoothL);
end
toc

csvwrite('F:\Clinical Gene expression network Project\Results\Data_smooth.txt', Data_smooth);   
save Data_smooth.mat
% SS=writetable( Data_smooth,'F:\Clinical Gene expression network Project\Results\Data_smooth.csv');


