
%% using kymograph_analysis data fit with Kelvin-Voigt
% y_data = a1(1-exp(-x_data/a2)
% Mimi suggested using lsqcurve fit for this
%       ahat: [2x1 double] WANT
%      resnorm: 1.2609e+03 WANT (SSE)
%     residual: [50x1 double] WANT
%     exitflag: 3
%       output: [1x1 struct]
%       lambda: [1x1 struct]
%     jacobian: [50x2 double]
for i = 1:3
    clear dispTemp temp
    if i == 1
        dispTemp = dispTS(:,tAblation:end);
    elseif i == 2
        dispTemp = dispTA(:,tAblation:end);
    else
        dispTemp = dispAE(:,tAblation:end);
    end
    x_input = timeMat(tAblation:end).';
    
    for j = 1:size(dispTemp,1)
        y_input = dispTemp(j,:).';
        s = size(x_input,1);
        % define the fitting function as an anon. function
        fit_KV = @(param,x_input) param(1)*(ones(s,1)-exp(-x_input./param(2)));
        param0 = [1;1]; % initial paramerter estimate
        % perform fitting with lsqcurvefit
        [temp(j).ahat,temp(j).resnorm,temp(j).residual,temp(j).exitflag,...
            temp(j).output,temp(j).lambda,temp(j).jacobian] = lsqcurvefit(fit_KV,param0,x_input,y_input);
        temp(j).TSS = sum((y_input-mean(y_input)).^2);
    end
    
    if i == 1
        ahat_TS = [temp.ahat];
        RSS_TS = [temp.resnorm];
        TSS_TS = [temp.TSS];
        residual_TS = [temp.residual];       
    elseif i == 2
        ahat_TA = [temp.ahat];
        RSS_TA = [temp.resnorm];
        TSS_TA = [temp.TSS];
        residual_TA = [temp.residual];
    else
        ahat_AE = [temp.ahat];
        RSS_AE = [temp.resnorm];
        TSS_AE = [temp.TSS];
        residual_AE = [temp.residual];       
    end
end
%% Calculate R2 for fit
% R2 = 1 - RSS/TSS
R2_TS = 1 - RSS_TS./TSS_TS;
R2_TA = 1 - RSS_TA./TSS_TA;
R2_AE = 1 - RSS_AE./TSS_AE;
%% Calculate v0
% v0 = derivative at time 0
v0_TS = ahat_TS(1,R2_TS > 0.5)./ahat_TS(2,R2_TS > 0.5);
v0_TA = ahat_TA(1,R2_TA > 0.5)./ahat_TA(2,R2_TA > 0.5);
v0_AE = ahat_AE(1,R2_AE > 0.5)./ahat_AE(2,R2_AE > 0.5);

%% Calculate time decay (tau)
tau_TS = ahat_TS(2,R2_TS > 0.5);
tau_TA = ahat_TA(2,R2_TA > 0.5);
tau_AE = ahat_AE(2,R2_AE > 0.5);
tau_AE(2) = [];

% calculate max displacment (a)
maxDisp_TS = ahat_TS(1,R2_TS > 0.5);
maxDisp_TA = ahat_TA(1,R2_TA > 0.5);
maxDisp_AE = ahat_AE(1,R2_AE > 0.5);
maxDisp_AE(2) = [];
%% Put into one matrix together and make boxplot of v0
v0_all = horzcat(v0_TS.', [v0_TA,NaN(1,4)].');
v0_all = horzcat(v0_all, [v0_AE,NaN(1,5)].');
figure
boxplot(v0_all,'labels',{'sqh-TS', 'sqh-TA', 'sqh-AE'})

%% wilcoxon rank sum test
p_ranksum = NaN(3,1);
p_ranksum(1) = ranksum(v0_TS, v0_TA);
p_ranksum(2) = ranksum(v0_TS, v0_AE);
p_ranksum(3) = ranksum(v0_TA, v0_AE);
%% ttest2
p_ttest2 = NaN(3,1);
[h, p_ttest2(1)] = ttest2(v0_TS, v0_TA);
[h, p_ttest2(2)] = ttest2(v0_TS, v0_AE);
[h, p_ttest2(3)] = ttest2(v0_TA, v0_AE);

%% put others into one matrix together and make boxplots
figure
tau_all = horzcat(tau_TS.', [tau_TA,NaN(1,4)].');
tau_all = horzcat(tau_all, [tau_AE,NaN(1,6)].');
maxDisp_all = horzcat(maxDisp_TS.', [maxDisp_TA,NaN(1,4)].');
maxDisp_all = horzcat(maxDisp_all, [maxDisp_AE,NaN(1,6)].');
% subplot(1,2,1)
% boxplot(tau_all,'labels',{'sqh-TS', 'sqh-TA', 'sqh-AE'})
% %h=findobj(gca,'tag','Outliers');
% %delete(h)
% ylim([0 30])
% subplot(1,2,2)
boxplot(maxDisp_all,'labels',{'sqh-TS', 'sqh-TA', 'sqh-AE'})
% ylim([0 9])
%% wilcoxon rank sum test for tau
p_ranksum_tau = NaN(3,1);
p_ranksum_tau(1) = ranksum(tau_TS, tau_TA);
p_ranksum_tau(2) = ranksum(tau_TS, tau_AE);
p_ranksum_tau(3) = ranksum(tau_TA, tau_AE);

p_ranksum_disp = NaN(3,1);
p_ranksum_disp(1) = ranksum(maxDisp_TS, maxDisp_TA);
p_ranksum_disp(2) = ranksum(maxDisp_TS, maxDisp_AE);
p_ranksum_disp(3) = ranksum(maxDisp_TA, maxDisp_AE);
%% ttest2
p_ttest2 = NaN(3,1);
[h, p_ttest2_tau(1)] = ttest2(tau_TS, tau_TA);
[h, p_ttest2_tau(2)] = ttest2(tau_TS, tau_AE);
[h, p_ttest2_tau(3)] = ttest2(tau_TA, tau_AE);
