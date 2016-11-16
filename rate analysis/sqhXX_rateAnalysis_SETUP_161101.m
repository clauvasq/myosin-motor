%% Load data from EDGE using load_edge_mat
% Need to make array that has folder with data, name of embryo, time
% resolution (seconds/frame), t-initial and t-final for segmentation.

% Initialize array with appropriate data:
% file name, XX, time res (sec/stack, start segmentation, end segmenation
sqhXXInfo={'1TS 20130213_2','TS',6.7,20,80;
    '1TS 20130213_4','TS',7,1,83;
    '1AE 20120821_1','AE',6.8,1,135;
    '1AE 20120822_2','AE',6.8,1,120;
    '1AE 20121031_1','AE',6.8,10,100;
    '1TA 20120816_3_CV','TA',6.7,25,150;
    '1TA 20120816_1','TA',6.7,1,73;
    '1TA 20121026_1','TA',6.7,25,90;
    '1AS 20120816_1','AS',7.5,45,150;
    '1AS 20121102_6', 'AS', 7.5, 1, 50;
    '1AS 20120816_2', 'AS', 7.5, 1, 43};

numEmbryos_sqhXX=size(sqhXXInfo,1);

%% Run load_edge_mat, adjusted for appropriate segmented times

sqhXXData=load_edge_mat(sqhXXInfo);

%% Calculate mean myosin, area, anisotropy and anistropy XY and SD for myo and area

sqhXXData=mean_area_myo_anisotropy(sqhXXData);

%% Plot mean myosin and area on same plot and then make subplot with all 
% areas from different embryos and myosins from different embryos on same
% plot

% Plot each embryo separately
% figure
% for i=1:numEmbryos_sqhXX
%     subplot(5,2,i)
%     plot_avg(sqhXXData(i).name,sqhXXData(i).meanArea,sqhXXData(i).meanMyo,...
%         sqhXXData(i).timeSec,'seconds');
% end

% Plot <myosin> for all embryos on same plot
cmp=colorcube(numEmbryos_sqhXX);
figure

hold on
for i=1:numEmbryos_sqhXX
    plot(sqhXXData(i).timeMin,...
        sqhXXData(i).meanArea,'Color',cmp(i,:));
end
% xlabel('time (sec)');
xlabel('time (min)');
ylabel('<area> (a.u.)');
legend(sqhXXData.name,'Location','SouthEastOutside');
hold off

% Plot <area> for all embryos on same plot
% figure
% subplot(1,2,2)
%%
figure
hold on
for i=1:numEmbryos_sqhXX
%     plot(sqhXXData(i).timeMin...
%         ,sqhXXData(i).meanArea, 'Color',cmp(i,:));
    plot(sqhXXData(i).meanArea);
end
% xlabel('time (sec)');
xlabel('time (min)');
ylabel('<Apical area> (um^2)');
legend(sqhXXData.name,'Location','SouthEastOutside');
hold off


%% Smooth data, normalize myosin to area, calculate rate and correlation 
% myosin and area

% choose time points to run cross-correlation calculation on
sqhXXData(1).corrTime=10:50;
sqhXXData(2).corrTime=20:60;
sqhXXData(3).corrTime=40:80;
sqhXXData(4).corrTime=80:120;
sqhXXData(5).corrTime=40:80;
sqhXXData(6).corrTime=10:50;
sqhXXData(7).corrTime=10:55;
sqhXXData(8).corrTime=20:60;
sqhXXData(9).corrTime=20:60;
sqhXXData(10).corrTime=1:40;
sqhXXData(11).corrTime=1:40;
% uses function smoo_rate_xcorr to smooth data, calculate rates and
% cross-correlation
sqhXXData=smoo_rate_xcorr(sqhXXData);

%% This script uses the rate of constriction to identify pulses, it 
% identifies the pulse as the maximum rate between two times that are
% greater than the specified threshold.

% load '~/Desktop/sqhXX 131220.mat'
for i=1:numEmbryos_sqhXX
    sqhXXData(i).numCells=size(sqhXXData(i).myo,2);
    sqhXXData(i).numFrames=size(sqhXXData(i).myo,1);
end

% Choose times of interest (300sec)
sqhXXPulse(1).endPts=[1 47];
sqhXXPulse(2).endPts=[20 64];
sqhXXPulse(3).endPts=[20 65];
sqhXXPulse(4).endPts=[30 75];
sqhXXPulse(5).endPts=[1 46];
sqhXXPulse(6).endPts=[10 56];
sqhXXPulse(7).endPts=[12 58];
sqhXXPulse(8).endPts=[5 51];
sqhXXPulse(9).endPts=[15 56];
sqhXXPulse(10).endPts=[1 41];
sqhXXPulse(11).endPts = [1 41];

figure
hold on
A = [1,3,7,10];
for i=9:11
% for j = 1:numEmbryos_sqhXX
%     i = A(j);
    t_I=sqhXXPulse(i).endPts(1);
    tF=sqhXXPulse(i).endPts(2);
    tFF=tF-t_I+1;
%     plot(sqhXXData(i).timeSec(1:tFF),sqhXXData(i).meanArea(t_I:tF));
    plot(sqhXXData(i).timeSec,sqhXXData(i).meanArea);
%     x = sqhXXData(i).timeSec(1:tFF);
%     y = sqhXXData(i).area(tI:tF);
%     shadedErrorBar(x, sqhXXData(i).meanArea(tI:tF), sqhXXData(i).stdArea(tI:tF));
end
xlabel('time (sec)');
ylabel('<Apical area> (um^2)');
% legend(sqhXXData.name,'Location','SouthEastOutside');
hold off

%% Choose values for threshold
thr=[0 0.5 1 2];

sqhXXPulse=pulse_id_thr(sqhXXData,sqhXXPulse,thr);

% min max normalize myosin
sqhXXPulse=myo_min_max_norm(sqhXXData,sqhXXPulse);
% sqhXXPulse=pulse_id_thrMyo(sqhXXData,sqhXXPulse,thr);
sqhXXPulse=pulse_id_thrAreaNormMyo(sqhXXData,sqhXXPulse,thr);

%% calculate time it takes to go from Ao to 0.5Ao
for i = 1:numEmbryos_sqhXX
    t_I = sqhXXPulse(i).endPts(1);
    area(i).initial = sqhXXData(i).area(t_I,:);
    t_end = sqhXXPulse(i).endPts(2);
%     area(i).pct_25 = repmat(0.25*area(i).initial,t_end,1);    
%     area(i).t_25 = area(i).pct_25 >= sqhXXData(i).area;
    area(i).pct_50 = repmat(0.50*area(i).initial,t_end-t_I+1,1);    
    area(i).t_50 = area(i).pct_50 >= sqhXXData(i).area(t_I:t_end,:);
    A = area(i).t_50;
    [I,J] = find(A==1);
    [~,m] = unique(J, 'first');
    area(i).index_50 = I(m);
    area(i).time_50 = (area(i).index_50 -tI)*sqhXXData(i).timeRes/60;
end
%%

sqhTS_time50 = vertcat(area(1).time_50,area(2).time_50);
sqhTA_time50 = vertcat(area(3).time_50,area(4).time_50, area(5).time_50);
sqhAE_time50 = vertcat(area(6).time_50,area(7).time_50, area(8).time_50);
sqhAS_time50 = vertcat(area(9).time_50,area(10).time_50, area(11).time_50);

%%
sqh_time50 = horzcat(sqhTS_time50,vertcat(sqhTA_time50, NaN(38,1)));
sqh_time50 = horzcat(sqh_time50, vertcat(sqhAE_time50, NaN(2,1)));
sqh_time50 = horzcat(sqh_time50, vertcat(sqhAS_time50, NaN(21,1)));
boxplot(sqh_time50,'labels',{'sqh-TS', 'sqh-TA', 'sqh-AE','sqh-AS'})
%%
p_ranksum_time50 = NaN(6,1);
p_ranksum_time50(1) = ranksum(sqhTS_time50, sqhTA_time50);
p_ranksum_time50(2) = ranksum(sqhTS_time50, sqhAS_time50);
p_ranksum_time50(3) = ranksum(sqhTS_time50, sqhAE_time50);
p_ranksum_time50(4) = ranksum(sqhTA_time50, sqhAS_time50);
p_ranksum_time50(5) = ranksum(sqhTA_time50, sqhAE_time50);
p_ranksum_time50(6) = ranksum(sqhAS_time50, sqhAE_time50);

%%
numCells = [sqhXXData(1).numCells+sqhXXData(2).numCells, ...
    sqhXXData(3).numCells + sqhXXData(4).numCells + sqhXXData(5).numCells,...
    sqhXXData(6).numCells + sqhXXData(7).numCells + sqhXXData(8).numCells,...
    sqhXXData(9).numCells + sqhXXData(10).numCells + sqhXXData(11).numCells];
pctCells = [size(sqhTS_time50,1), size(sqhTA_time50,1),...
    size(sqhAE_time50,1), size(sqhAS_time50,1)];
pctCells = pctCells./numCells;