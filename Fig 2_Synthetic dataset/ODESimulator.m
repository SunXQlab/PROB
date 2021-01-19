
% %%% Simulated Coeficient Matrix
% C=csvread('F:\Clinical Gene expression network Project\Code\Simulated ODE data - 6GeneNet\Coefficient Matrix_new.csv');
% C=csvread('F:\Clinical Gene expression network Project\Code\Simulated ODE data - 6GeneNet\Coefficient Matrix_double3net.csv');
C=csvread('F:\E Disk\Clinical Gene expression network Project\Code\Simulated ODE data - 6GeneNet\Coefficient Matrix_feedback.csv');

D=diag(C)';
C(logical(eye(size(C,1))))=0;
% D=D*0.8;
% C=C*0.5;

TimeLength=100;
X0=ones(1,size(C,1)); %.*rand(1,size(C,1));
% X0=[1 6 1];
X=zeros(size(C,1),TimeLength);
X(:,1)=X0;

h=0.01;
%     for j=1:TimeLength-1
%         X(:,j+1)=X(:,j)+C(:,:)*X(:,j)*h;  %  .*X(:,j)+D(:).*X(:,j)).*0.01
% %         X(X>1)=1; X(X<0)=0;
%     end
%     
% %%% Another method for solving ODEs    
  for tt=1:TimeLength-1
       for i=1:size(C,1)
             f=0;
            for j=1:size(C,1)
            f=f+C(i,j)*X(j,tt)*X(i,tt);
            end
             X(i,tt+1)=X(i,tt)+(f)*h;
        end
 end

    for j=1:TimeLength-1
        X(:,j+1)=X(:,j)+((C(:,:)*X(:,j)).*X(:,j)+D(:).*X(:,j)).*0.01;  %
%         X(X>1)=1; X(X<0)=0;
    end

    %%% Another method for solving ODEs    
  for tt=1:TimeLength-1
       for i=1:size(C,1)
             f=0;
            for j=1:size(C,1)
            f=f+C(i,j)*X(j,tt)*X(i,tt);
            end
             X(i,tt+1)=X(i,tt)+(f+D(i)*X(i,tt))*h;
        end
  end

% X=X./(max(X,[],2))
% figure,
% pcolor(X')

figure,
hold on, plot(h:h:1,X(1,:),'-','color',[0.5 0.5 0.5],'LineWidth',3)
hold on, plot(h:h:1,X(2,:),'-','color',[0.5 0.5 1],'LineWidth',3)
hold on, plot(h:h:1,X(3,:),'-','color',[0.5 1 0.5],'LineWidth',3)
hold on, plot(h:h:1,X(4,:),'-','color',[1 0.5 0.5],'LineWidth',3)
hold on, plot(h:h:1,X(5,:),'-','color',[0.8 0.3 0.5],'LineWidth',3)
hold on, plot(h:h:1,X(6,:),'-','color',[0.5 0.3 0.8],'LineWidth',3)
% hold on, plot(h:h:1,X(7,:),'-','color',[0.3 0.3 1],'LineWidth',3)
% hold on, plot(h:h:1,X(8,:),'-','color',[0.3 1 0.8],'LineWidth',3)
% hold on, plot(h:h:1,X(9,:),'-','color',[0.7 0.5 0.8],'LineWidth',3)
% hold on, plot(h:h:1,X(10,:),'-','color',[0.6 0.5 0.3],'LineWidth',3)
set(gca,'FontSize',15)
xlabel('Time','FontSize',20);
ylabel('Gene expression','FontSize',20)
legendlabel=legend('G1','G2','G3','G4','G5','G6','Location','NorthEastoutside','Box','off')  %'G7','G8','G9','G10',
set(legendlabel,'Box','off')

%%%%%% Simulate clinical gene expression data
rng('shuffle') 

%% Guassian noise
% NoiseLevel=10;
% X_noise=X.*(1.+NoiseLevel*1e-2.*randn(size(X)));  % 5% noise
%% drop-out noise
% DropoutNoise=randsample([0 1],600,true,[0.05 0.95]);
% DropoutNoise=reshape(DropoutNoise,6,100);
% X_noise=X.*DropoutNoise;     % drop-out noise
%% Exponential noise
NoiseLevel=0.1;
E_Noise=exprnd(NoiseLevel,size(X));
X_noise=X.*(1.+E_Noise);
%%
Stage=[ones(1,25),ones(1,25)*2,ones(1,25)*3,ones(1,25)*4];
X_stage=[X_noise;Stage;1:100];
rand_index=randperm(100);
X_stage=X_stage(:,rand_index);
% save X_stage.mat
% xlswrite('F:\Clinical Gene expression network Project\Code\Simulated ODE data\GeneExpression_Stage_Simulated.xls', X_stage);  

figure,
hold on, plot(X_noise(1,:),'-','color',[0.5 0.5 0.5],'LineWidth',3)
hold on, plot(X_noise(2,:),'-','color',[0.5 0.5 1],'LineWidth',3)
hold on, plot(X_noise(3,:),'-','color',[0.5 1 0.5],'LineWidth',3)
hold on, plot(X_noise(4,:),'-','color',[1 0.5 0.5],'LineWidth',3)
hold on, plot(X_noise(5,:),'-','color',[0.8 0.3 0.5],'LineWidth',3)
hold on, plot(X_noise(6,:),'-','color',[0.8 0.3 0.8],'LineWidth',3)
set(gca,'FontSize',15)
xlabel('Real progression','FontSize',20);
ylabel('Gene expression','FontSize',20)
Legendlabel=legend('G1','G2','G3','G4','G5','G6','Location','NorthEastoutside','Box','off');  %'G7','G8','G9','G10',


figure,
hold on, plot(X_stage(1,:),'o','color',[0.5 0.5 0.5],'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5])
hold on, plot(X_stage(2,:),'o','color',[0.5 0.5 1],'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',[0.5 0.5 1])
hold on, plot(X_stage(3,:),'o','color',[0.5 1 0.5],'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',[0.5 1 0.5])
hold on, plot(X_stage(4,:),'o','color',[1 0.5 0.5],'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',[1 0.5 0.5])
hold on, plot(X_stage(5,:),'o','color',[0.8 0.3 0.5],'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',[0.8 0.3 0.5])
hold on, plot(X_stage(6,:),'o','color',[0.8 0.3 0.8],'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',[0.8 0.3 0.8])
hold on, plot(1:25,zeros(1,25)+min(min(X_stage(1:3,:)))-1,'-','color',[0.75 0.75 0.75],'LineWidth',8);hold on, plot(26:50,zeros(1,25)+min(min(X_stage(1:3,:)))-1,'-','color',[0.65 0.65 0.65],'LineWidth',8);;hold on, plot(51:75,zeros(1,25)+min(min(X_stage(1:3,:)))-1,'-','color',[0.5 0.5 0.5],'LineWidth',8);hold on, plot(76:100,zeros(1,25)+min(min(X_stage(1:3,:)))-1,'-','color',[0.4 0.4 0.4],'LineWidth',8);
hold on, plot(1:25,zeros(1,25)-2,'-','color',[0.75 0.75 0.75],'LineWidth',8);hold on, plot(26:50,zeros(1,25)-2,'-','color',[0.65 0.65 0.65],'LineWidth',8);hold on, plot(51:75,zeros(1,25)-2,'-','color',[0.5 0.5 0.5],'LineWidth',8);hold on, plot(76:100,zeros(1,25)-2,'-','color',[0.4 0.4 0.4],'LineWidth',8);
set(gca,'FontSize',15)
xlabel('Randomized sample ID','FontSize',20);
ylabel('Gene expression','FontSize',20)
set(gca,'xtick',[12.5:25:100]);
set(gca,'xticklabel',{'S1','S2','S3','S4'});
Leg=legend('G1','G2','G3','G4','G5','G6','Location','NorthEastOutside') %,'G4','G5','Location','NorthWest','Box','off')
set(Leg,'Box','off')
% axis([0 100 min(min(X_stage(1:3,:)))-0.4 max(max(X_stage(1:3,:)))+0.4])
axis([0 100 -1 8])


% %%%% 
% X_heatmap=[X_stage(1:6,X_stage(7,:)==1),X_stage(1:6,X_stage(7,:)==2),X_stage(1:6,X_stage(7,:)==3),X_stage(1:6,X_stage(7,:)==4)];
% figure, heatmap(X_heatmap);
% % X_heatmap=X_stage;
% Size=30;
% Color=X_heatmap(4,:);
% ExpGene=X_heatmap(1:3,:);
% figure,scatter3(ExpGene(1,:),ExpGene(2,:),ExpGene(3,:),Size,Color,'filled')
% xlabel('Gene A','Rotation',15)
% ylabel('Gene B','Rotation',-18)
% zlabel('Gene C')
% %%% color: 1: dark blue; 2: shallow blue; 3: yellow; 4: green.

