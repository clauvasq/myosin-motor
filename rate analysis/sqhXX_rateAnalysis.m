% instants of constriction rates with comparable myosin decreases

%% Calculate cutoff using TS embryos and then use for all embryos.
thr = [1,2];
smooRateMyoTS=horzcat(sqhXXPulse(1).smooRateMyo,vertcat(sqhXXPulse(2).smooRateMyo,NaN(2,67)));
meanRateMyoTS=nanmean(smooRateMyoTS(:));
stdRateMyoTS=nanstd(smooRateMyoTS(:));
cutoff=thr*stdRateMyoTS+meanRateMyoTS;

%% Find where smoothed rate Myo for each embryo is greater than 1 s.d. 
% above the mean, but less than 2 s.d.
for k=1:numEmbryos_sqhXX
    A = sqhXXPulse(k).smooRateMyo;
    L = cutoff(1);
    U = cutoff(2);
    sqhXXPulse(k).bound = A(A>L & A<U);    
    sqhXXPulse(k).boundArea = sqhXXPulse(k).smooRateArea(A>L & A<U);
    sqhXXPulse(k).numBound = size(sqhXXPulse(k).bound, 1);
    sqhXXPulse(k).thr2 = A(A>U);
    sqhXXPulse(k).thr2Area = sqhXXPulse(k).smooRateArea(A>U);
    sqhXXPulse(k).numThr2 = size(sqhXXPulse(k).thr2, 1);
end

%% Plot histograms
figure
for i = 1:numEmbryos_sqhXX
%     subplot(2,5,i)
    subplot(4,3,i)
    rateArea(i).histProp = histogram(sqhXXPulse(i).boundArea,...
        'Normalization','probability');
    title(sqhXXData(i).name);
end

%% Pool mutants together and plot
sqhTSRateArea = vertcat(sqhXXPulse(1).boundArea, sqhXXPulse(2).boundArea);
sqhAERateArea = vertcat(sqhXXPulse(3).boundArea, sqhXXPulse(4).boundArea, sqhXXPulse(5).boundArea);
sqhTARateArea = vertcat(sqhXXPulse(6).boundArea, sqhXXPulse(7).boundArea, sqhXXPulse(8).boundArea);
sqhASRateArea = vertcat(sqhXXPulse(9).boundArea, sqhXXPulse(10).boundArea, sqhXXPulse(11).boundArea);
% sqhASRateArea = vertcat(sqhXXPulse(10).boundArea);


% sqhTSRateArea_thr2 = vertcat(sqhXXPulse(1).thr2Area, sqhXXPulse(2).thr2Area);
% sqhAERateArea_thr2 = vertcat(sqhXXPulse(3).thr2Area, sqhXXPulse(4).thr2Area, sqhXXPulse(5).thr2Area);
% sqhTARateArea_thr2 = vertcat(sqhXXPulse(6).thr2Area, sqhXXPulse(7).thr2Area, sqhXXPulse(8).thr2Area);
% sqhASRateArea_thr2 = vertcat(sqhXXPulse(9).thr2Area, sqhXXPulse(10).thr2Area);


figure
hold on
histogram(sqhTSRateArea, 'Normalization','probability');
histogram(sqhTARateArea, 'Normalization','probability');
histogram(sqhAERateArea, 'Normalization','probability');
histogram(sqhASRateArea, 'Normalization','probability');
legend( {'TS', 'TA', 'AE', 'AS'});
hold off

%% Make boxplot
% uses area normalized myosin
allRateArea = vertcat(sqhTSRateArea, NaN(292,1));
allRateArea = horzcat(allRateArea, sqhTARateArea);
tempRateArea = vertcat(sqhAERateArea, NaN(365, 1));
allRateArea = horzcat(allRateArea, tempRateArea);
tempRateArea = vertcat(sqhASRateArea, NaN(319, 1));
allRateArea = horzcat(allRateArea, tempRateArea);

% %
% allRateArea = vertcat(sqhTSRateArea, NaN(326,1));
% allRateArea = horzcat(allRateArea, sqhTARateArea);
% tempRateArea = vertcat(sqhAERateArea, NaN(261, 1));
% allRateArea = horzcat(allRateArea, tempRateArea);
% tempRateArea = vertcat(sqhASRateArea, NaN(331, 1));
% allRateArea = horzcat(allRateArea, tempRateArea);

% uses thr2, arean normalized
% allRateArea_thr2 = vertcat(sqhTSRateArea_thr2, NaN(188,1));
% allRateArea_thr2 = horzcat(allRateArea_thr2, sqhTARateArea_thr2);
% tempRateArea = vertcat(sqhAERateArea_thr2, NaN(357, 1));
% allRateArea_thr2 = horzcat(allRateArea_thr2, tempRateArea);
% tempRateArea = vertcat(sqhASRateArea_thr2, NaN(261, 1));
% allRateArea_thr2 = horzcat(allRateArea_thr2, tempRateArea);

% % uses thr2
% allRateArea_thr2 = vertcat(sqhTSRateArea_thr2, NaN(144,1));
% allRateArea_thr2 = horzcat(allRateArea_thr2, sqhTARateArea_thr2);
% tempRateArea = vertcat(sqhAERateArea_thr2, NaN(250, 1));
% allRateArea_thr2 = horzcat(allRateArea_thr2, tempRateArea);
% tempRateArea = vertcat(sqhASRateArea_thr2, NaN(162, 1));
% allRateArea_thr2 = horzcat(allRateArea_thr2, tempRateArea);


figure
boxplot(allRateArea, 'labels', {'sqh-TS', 'sqh-TA', 'sqh-AE','sqh-AS'});
h=findobj(gca,'tag','Outliers');
% set(h,'Visible', 'off')

 %% and do ranksum test
 pNorm = NaN(4,1);
[h pNorm(1)] = kstest(sqhTSRateArea);
[h pNorm(2)] = kstest(sqhTARateArea);
[h pNorm(3)] = kstest(sqhAERateArea);
[h pNorm(4)] = kstest(sqhASRateArea);

 p_ranksum = NaN(4,1);
 p_ranksum(1) = ranksum(sqhTSRateArea, sqhTARateArea);
 p_ranksum(2) = ranksum(sqhTSRateArea, sqhAERateArea);
 p_ranksum(3) = ranksum(sqhTSRateArea, sqhASRateArea);
 p_ranksum(4) = ranksum(sqhTARateArea, sqhAERateArea);
 
p_ttest2 = NaN(4,1);
[h, p_ttest2(1)] = ttest2(sqhTSRateArea, sqhTARateArea);
[h, p_ttest2(2)] = ttest2(sqhTSRateArea, sqhAERateArea);
[h, p_ttest2(3)] = ttest2(sqhTSRateArea, sqhASRateArea);
[h, p_ttest2(4)] = ttest2(sqhTARateArea, sqhAERateArea);

%% cdfplot
figure
hold on
cdfplot(sqhTSRateArea);
cdfplot(sqhTARateArea);
cdfplot(sqhAERateArea);
cdfplot(sqhASRateArea);
legend({'TS', 'TA', 'AE', 'AS'});
hold off
 