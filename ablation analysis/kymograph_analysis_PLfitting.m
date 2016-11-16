


 %% Fit power law to  displacement curves 
%  y = K[x^a] , where b = log(K)
% log(y) = alog(x) + b
% a = .PL(j).coeff(1) [THIS IS THE VALUE I'M INTERESTED IN]
% b = .PL(j).coeff(2)

delta = 1;
deltaEnd = 33;
for i = 1:3
% for i =1
    clear dispTemp tempPL
    if i == 1
        dispTemp = dispTS(:,tAblation+delta:end-deltaEnd);
    elseif i == 2
        dispTemp = dispTA(:,tAblation+delta:end-deltaEnd);
    else
        dispTemp = dispAE(:,tAblation+delta:end-deltaEnd);
    end
    x_input = timeMat(tAblation+delta:end-deltaEnd);
%     x_input(1) = 0.01;

    
    for j = 1:size(dispTemp,1)
%     for j = 1
        y_input = dispTemp(j,:);
        [tempPL(j).coeff,tempPL(j).s] = polyfit(log(x_input), log(y_input),1);
        tempPL(j).R2 = R2_polyfit(log(x_input),log(y_input), tempPL(j).coeff);
    end
    
    if i == 1
        PLcoeff_TS = [tempPL.coeff]; 
        PL_R2_TS = [tempPL.R2];
    elseif i == 2
        PLcoeff_TA = [tempPL.coeff]; 
        PL_R2_TA = [tempPL.R2];
    else
        PLcoeff_AE = [tempPL.coeff];
        PL_R2_AE = [tempPL.R2];
    end
end

%% get alpha value (odd values of PLcoeff_XX matrices)
alpha_TS = real(PLcoeff_TS(1:2:end));
% alpha_TS = alpha_TS(PL_R2_TS > 0.5);
alpha_TA = real(PLcoeff_TA(1:2:end));
% alpha_TA = alpha_TA(PL_R2_TA > 0.5);
alpha_AE = real(PLcoeff_AE(1:2:end));
% alpha_AE = alpha_AE(PL_R2_AE > 0.5);

beta_TS = real(PLcoeff_TS(2:2:end));
beta_TA = real(PLcoeff_TA(2:2:end));
beta_AE = real(PLcoeff_AE(2:2:end));

%%
figure
subplot(1,3,1)

hold on
plot(log(x_input), log(dispTS(:, tAblation+delta:end-deltaEnd)).','o')
for i = 1:n
    plot(log(x_input), alpha_TS(i)*log(x_input)+beta_TS(i))
end
title('sqh-TS')

subplot(1,3,2)
hold on
plot(log(x_input), log(dispTA(:, tAblation+delta:end-deltaEnd)).','o')
for i = 1:size(alpha_TA,2)
    plot(log(x_input), alpha_TA(i)*log(x_input)+beta_TA(i))
end
title('sqh-TA')

subplot(1,3,3)
hold on
plot(log(x_input), log(dispAE(:, tAblation+delta:end-deltaEnd)).','o')
for i = 1:size(alpha_AE,2)
    plot(log(x_input), alpha_AE(i)*log(x_input)+beta_AE(i))
end
title ('sqh-AE')

%% get rid of values that are less than zero or greater than 1
beta_AE = beta_AE(alpha_AE>=0);
alpha_AE = alpha_AE(alpha_AE>=0);
beta_TA = beta_TA(alpha_TA<=1);
alpha_TA = alpha_TA(alpha_TA<=1);
%% make boxplot of alphas
alpha_all = horzcat(alpha_TS.', [alpha_TA, NaN(1,5)].');
alpha_all = horzcat(alpha_all, [alpha_AE, NaN(1,7)].');
figure
boxplot(alpha_all,'labels',{'sqh-TS', 'sqh-TA', 'sqh-AE'})
ylabel('alpha')
h = findobj(gcf,'tag','Outliers');
xdata = get(h,'XData');
ydata = get(h,'YData');

%%
pAlpha_ranksum = NaN(3,1);
pAlpha_ranksum(1) = ranksum(alpha_TS, alpha_TA);
pAlpha_ranksum(2) = ranksum(alpha_TS, alpha_AE);
pAlpha_ranksum(3) = ranksum(alpha_TA, alpha_AE);

