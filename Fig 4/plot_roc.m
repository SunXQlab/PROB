     function  ROC_data = plot_roc( predict, ground_truth, dispp, dispt )  
    % ------------------------------------------------------------------- %
    % INPUTS  
    %  predict       - 分类器对测试集的分类结果  
    %  ground_truth - 测试集的正确标签,这里只考虑二分类，即0和1  
% ----------------------------------------------------------------------- %                            
    % OUTPUTS      
%       * dispp:    (Optional) If dispp is 1, the ROC Curve will be disp- %
%                   ayed inside the active figure. If dispp is 0, no figure
%                   will be displayed.                                    %
%       * dispt:    (Optional) If dispt is 1, the optimum threshold para- %
%                   meters obtained will be displayed on the MATLAB log.  %
%                   Otherwise, if dispt is 0, no parameters will be disp- %
%                   ayed there.                                           %
                                       %
% ----------------------------------------------------------------------- %                                                                       %
%   Output variables                                                      %
%       * ROC_data: Struct that contains all the curve parameters.        %
%           - param:    Struct that contains the cuantitative parameters  %
%                       of the obtained curve, which are:                 %
%               + Threshold:Optimum threshold calculated in order to maxi-%
%                           mice the sensitivity and specificity values,  %
%                           which is colocated in the nearest point to    %
%                           (0,1).                                        %
%               + AROC:     Area under ROC curve.                         %
%               + Accuracy: Accuracy.                                     %
%               + Sensi:    Sensitivity (i.e., recall, hit rate, or true  %
%                           positive rate).                               %
%               + Speci:    Specificity (i.e., selectivity, or true       %
%                           negative rate).                               %
%               + PPV:      Positive predicted value (i.e., precision).   %
%               + NPV:      Negative predicted value.                     %
%               + FNR:      False negative rate (i.e., miss rate).        %
%               + FPR:      False positive rate (i.e., fall-out).         %
%               + FDR:      False discovery rate.                         %
%               + FOR:      False omission rate.                          %
%               + F1_score: F1 score (harmonic mean of precision and      %
%                           sensitivity).                                 %
%               + MCC:      Matthews correlation coefficient.             %
%               + BM:       Bookmaker informedness.                       %
%               + MK:       Markedness.                                   %
%           - curve:    Matrix that contains the specificity and specifi- %
%                       city of each threshold point in columns.       
% ----------------------------------------------------------------------- %
%   Author: Xiaoqiang Sun @ SYSU                                      %
%   Date:    8/1/2019                                                   %
%   E-mail:  xiaoqiangsun88 at gmail dot com  
% ----------------------------------------------------------------------- %

    C=ground_truth;
    S=predict;
    
    predict=reshape(predict,1,size(predict,1)*size(predict,2));  
    ground_truth=reshape(ground_truth,1,size(ground_truth,1)*size(ground_truth,2));  
    
    curve = zeros(length(ground_truth),2);  % Sensitivity
    distance = zeros(length(ground_truth),1);  % Specificity
    
    %初始点为（1.0, 1.0）  
    x = 1.0;  
    y = 1.0;  
    %计算出ground_truth中正样本的数目pos_num和负样本的数目neg_num  
    pos_num = sum(ground_truth==1);  
    neg_num = sum(ground_truth==0);  
    %根据该数目可以计算出沿x轴或者y轴的步长  
    x_step = 1.0/neg_num;  
    y_step = 1.0/pos_num;  
    %首先对predict中的分类器输出值按照从小到大排列  
    [predict,index] = sort(predict);  
    ground_truth = ground_truth(index);  
    %对predict中的每个样本分别判断他们是FP或者是TP  
    %遍历ground_truth的元素，  
    %若ground_truth[i]=1,则TP减少了1，往y轴方向下降y_step  
    %若ground_truth[i]=0,则FP减少了1，往x轴方向下降x_step  
    for i=1:length(ground_truth)  % Threshold
        if ground_truth(i) == 1  
            y = y - y_step;  
        else  
            x = x - x_step;  
        end  
        TP(i)=y;  
        FP(i)=x;  
        
    end  
    
%     %%%%% TN  FN 
%     x=0;
%     y=0;
%     x_step = 1.0/pos_num;  
%     y_step = 1.0/ neg_num;  
%     [predict,index] = sort(predict);  
%     ground_truth = ground_truth(index);  
%     %对predict中的每个样本分别判断他们是FP或者是TP  
%     %遍历ground_truth的元素，  
%     %若ground_truth[i]=1,则TP减少了1，往y轴方向下降y_step  
%     %若ground_truth[i]=0,则FP减少了1，往x轴方向下降x_step  
%     for i=1:length(ground_truth)  
%         if ground_truth(i) == 1  
%             x= x + x_step;    
%         else  
%             y = y + y_step; 
%         end  
%         TN(i)=x;  
%         FN(i)=y;  
%     end  
    
    %%% More simply...
    FN=1-TP;
    TN=1-FP;
    
      % Calculating the sensibility and specificity of each threshold
        curve(:,1) = 1-TN./(TN + FP);  % 1-Specificity 
        curve(:,2) = TP./(TP + FN);	% Sensitivity
        
        % Distance between each point and the optimum point (0,1)
        distance(:)= sqrt((curve(:,1)).^2+(curve(:,2)-1).^2);
        
        % Optimum threshold and parameters
    [ind, opt] = min(distance);
    TP = TP(opt);    % No. true positives
    FP = FP(opt);    % No. false positives 
    FN = FN(opt);     % No. false negatives                                 
    TN = TN(opt);     % No. true negatives      
      
    
    %  calculate AUC
    labels=reshape(C~=0,size(C,1)*size(C,2),1);
    scores=reshape(S,size(S,1)*size(S,2),1);
    AUC =scoreAUC(labels, scores);
    param.AROC = AUC;
    
    % Output parameters
    param.Threshold = predict(opt);               % Optimum threshold position
    param.Sensi = curve(opt,2);                 % Sensitivity
    param.Speci = 1-curve(opt,1);                 % Specificity
    param.AROC  = AUC; % Area under curve
    param.Accuracy = (TP+TN)/(TP+TN+FP+FN);     % Aaccuracy
    param.PPV   = TP/(TP+FP);                   % Positive predictive value
    param.NPV   = TN/(TN+FN);                   % Negative predictive value
    param.FNR   = FN/(FN+TP);                   % False negative rate
    param.FPR   = FP/(FP+TN);                   % False positive rate
    param.FDR   = FP/(FP+TP);                   % False discovery rate
    param.FOR   = FN/(FN+TN);                   % False omission rate
    param.F1_score = 2*TP/(2*TP+FP+FN);         % F1 score
    param.F2_score = 5*TP/(4*TP+FP+FN);         % F1 score
    param.MCC   = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));  % Matthews correlation coefficient
    param.BM    = param.Sensi+param.Speci-1;    % Informedness
    param.MK    = param.PPV+param.NPV-1;        % Markedness
    
    param.TP = TP;    % No. true positives
    param.FP = FP;    % No. false positives 
    param.FN = FN;    % No. false negatives                                 
    param.TN = TN;    % No. true negatives  
    
    %画出图像 
%     figure,
%     plot(TP,FP,'-ro','LineWidth',2,'MarkerSize',3);  
%     xlabel('FPR');  
%     ylabel('TPR');  
%     title('ROC curve');
    
      if(dispp == 1)
%         figure,
        fill_color = [11/255, 208/255, 217/255];
        curve_new=[[1,1];curve];
%         fill(curve(:,1),curve(:,2),  fill_color,'FaceAlpha',0.5);
  
        hold on; plot(curve_new(:,1),curve_new(:,2),  '.-','color',[rand,rand,rand], 'LineWidth', 2,'MarkerSize',2);  % 
%         hold on; plot(curve(opt,1),curve(opt,2),  'or', 'MarkerSize', 10);
%         hold on; plot(curve(opt,1),curve(opt,2),  'xr', 'MarkerSize', 12);
%         hold on; plot(0:0.01:1,0:0.01:1,  '-.','color',[0.50,0.54,0.53], 'LineWidth', 1);
        hold off; axis square; grid on; xlabel('False positive rate'); ylabel('True positive rate');
        set(gca,'FontSize',15)
        axis([0 1 0 1])
        title(['AROC = ' num2str(param.AROC)]);
      end
    
    
    % Log screen parameters if required
    if(dispt == 1)
        fprintf('\n ROC CURVE PARAMETERS\n');
        fprintf(' ------------------------------\n');
        fprintf('  - Distance:     %.4f\n', distance(opt));
        fprintf('  - Threshold:    %.4f\n', param.Threshold);
        fprintf('  - Sensitivity:  %.4f\n', param.Sensi);
        fprintf('  - Specificity:  %.4f\n', param.Speci);
        fprintf('  - AROC:         %.4f\n', param.AROC);
        fprintf('  - Accuracy:     %.4f\n', param.Accuracy);
        fprintf('  - PPV:          %.4f\n', param.PPV);
        fprintf('  - NPV:          %.4f\n', param.NPV);
        fprintf('  - FNR:          %.4f\n', param.FNR);
        fprintf('  - FPR:          %.4f\n', param.FPR);
        fprintf('  - FDR:          %.4f\n', param.FDR);
        fprintf('  - FOR:          %.4f\n', param.FOR);
        fprintf('  - F1 score:     %.4f\n', param.F1_score);
        fprintf('  - F2 score:     %.4f\n', param.F2_score);
        fprintf('  - MCC:          %.4f\n', param.MCC);
        fprintf('  - BM:           %.4f\n', param.BM);
        fprintf('  - MK:           %.4f\n', param.MK);
        fprintf(' \n');
    end
    
    % Assinging parameters and curve data
    ROC_data.param = param;
    ROC_data.curve = curve;
    end  