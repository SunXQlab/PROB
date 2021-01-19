%% reform the score matrix to be compared with [D',C] matrix
clear S_new
% clear C_new
% S=(S-min(S,[],2))./(max(S,[],2)-min(S,[],2));
S_new=S;
for i=1:size(C,1)
    S_new(i,i)=0;
    S_new(i,setdiff(1:size(X,1),i))=S(i,1:end-1);
end
C(logical(eye(size(C,1))))=0;

% S_new=[];
% C_new=[];
% for i=1:size(X,1)
%     S_new=[S_new,S(i,setdiff(1:size(X,1),i))];
%     C_new=[C_new,C(i,setdiff(1:size(X,1),i))];
% end


% ROC_data=plot_roc(S_new, [C,D']~=0 ,1,1);  % [C,D'] correspond to Intercept term and C 
% ROC_data=plot_roc(S_new, C_new~=0 ,1,1);  % [C,D'] correspond to Intercept term and C 
ROC_data=plot_roc(S_new, C~=0 ,1,1);  % [C,D'] correspond to Intercept term and C 

% S_normalized=(S-min(S,[],2))./(max(S,[],2)-min(S,[],2))
% S_normalized=zscore(S')';
% plot_roc( S_normalized, C~=0 ,1, 1)  % 1; 
% 
% 
% 
% labels=reshape(C~=0,size(X,1)*size(X,1),1);
% scores=reshape(S,size(S,1)*size(S,2),1);
% 
% AUC =scoreAUC(labels, scores);
% 
% %  [~,~,~,AUC2] = perfcurve(labels,scores,posclass); %  Another method for calculating AUAC


for i=1:size(X,1)
    CC(i,size(X,1))=D(i);
    CC(i,setdiff(1:size(X,1),size(X,1)))=C(i,setdiff(1:size(X,1),i));
end

Para_new=zeros(size(X,1),size(X,1),1e4);
for i=1:size(C,1)
    Para_new(i,i,:)=zeros(1,1,1e4);
    Para_new(i,setdiff(1:size(X,1),i),:)=Para_Post_pdf{i}.BetaDraws(1:end-1,:);
end

figure,
fig_ind=1;
for i=1:size(X,1)
    for j=1:size(X,1)
        subplot(size(X,1),size(X,1),fig_ind)
        [f,xi]=ksdensity(reshape(Para_new(i,j,:),1,1e4));  
        set(gca,'LineWidth',1)
        ksdensity(reshape(Para_new(i,j,:),1,1e4)); 
        hold on; plot(C(i,j)*ones(1,length(0:1e-2:max(f)+0.2)),0:1e-2:max(f)+0.2,'r-','LineWidth',1);
        set(gca,'FontSize',10)
        axis([min(C(i,j),min(xi)) max(C(i,j),max(xi)) 0 max(f)+0.2])
        fig_ind=fig_ind+1;
    end
end

