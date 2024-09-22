%读取数据
% M = readtable('C:\Users\Administrator\Desktop\华为杯竞赛\datas\附件一（训练集）.xlsx');
% 
% figue(1);
% signal = table2array(M(2, 5:end));
% plot(signal);
% 
% signal2 = table2array(M(3, 5:end));
% plot(signal);

clear;
clc;


%%读取文件夹下所有tif文件
%覆被数据
dirOutput = dir('C:/Users/Administrator/Desktop/华为杯竞赛/中国大陆0.5°土地利用和覆盖变化数据集(1900-2019年)/数据实体/*.tif');
% dirOutput = dir('C:/Users/Administrator/Desktop/华为杯竞赛/中国大陆0.5°土地利用和覆盖变化数据集(1900-2019年)/数据实体/cropland*.tif');
% dirOutput = dir('C:/Users/Administrator/Desktop/华为杯竞赛/中国大陆0.5°土地利用和覆盖变化数据集(1900-2019年)/数据实体/forest*.tif');
% dirOutput = dir('C:/Users/Administrator/Desktop/华为杯竞赛/中国大陆0.5°土地利用和覆盖变化数据集(1900-2019年)/数据实体/grass*.tif');
% dirOutput = dir('C:/Users/Administrator/Desktop/华为杯竞赛/中国大陆0.5°土地利用和覆盖变化数据集(1900-2019年)/数据实体/shrub*.tif');
% dirOutput = dir('C:/Users/Administrator/Desktop/华为杯竞赛/中国大陆0.5°土地利用和覆盖变化数据集(1900-2019年)/数据实体/wetland*.tif');

folder = string({dirOutput.folder}');
file = string({dirOutput.name}');
filepath = strcat(folder, '\', file);

%降水数据
%文件名
%读取NetC数据
ncdisp('CHM_PRE_0.25dg_19612022.nc');
rainData = ncread('CHM_PRE_0.25dg_19612022.nc', 'pre');
rainTime = ncread('CHM_PRE_0.25dg_19612022.nc', 'time');
rainYear = ncread('CHM_PRE_0.25dg_19612022.nc', 'years');
%% 问一处理覆被数据集
numFilePath = length(filepath);
yearsNum = 2020 - 1990;

totalM = ones(78,127);
figure(3);
nonZeroSum = zeros(5,yearsNum);%120为1990年到2020年的年份数

cropland_idx = 1;
forest_idx = 1;
grass_idx = 1;
shrub_idx = 1;
wetland_idx = 1;

for i = 1 : numFilePath 
    filepath_t = filepath{i};
    % disp(filepath_t);
    pause(0.0001);
    
    
    % if mod(i,yearsNum) == 0
    %     totalM = totalM - ansMatrix;
    % end
    
    
    if contains(file{i}, 'cropland')  
        yearStr = extractBetween(file{i}, 'cropland-', '.tif');%提取所需要的年份
        year = str2double(yearStr);
        if year >= 1990 && year <= 2020
            ansMatrix = readTIF(filepath_t, '农业用地', year);
            nonZeroSum(1,cropland_idx) = sum(ansMatrix(ansMatrix ~= 0));
            cropland_idx = cropland_idx + 1;
        end

    elseif contains(file{i}, 'forest')
        yearStr = extractBetween(file{i}, 'forest-', '.tif');%提取所需要的年份
        year = str2double(yearStr);
        if year >= 1990 && year <= 2020
            ansMatrix = readTIF(filepath_t, '森林', year);
            nonZeroSum(2,forest_idx) = sum(ansMatrix(ansMatrix ~= 0));
            forest_idx = forest_idx + 1;
        end
    elseif contains(file{i}, 'grass')
        yearStr = extractBetween(file{i}, 'grass-', '.tif');%提取所需要的年份
        year = str2double(yearStr);
        if year >= 1990 && year <= 2020
            ansMatrix = readTIF(filepath_t, '草原', year);
            nonZeroSum(3,grass_idx) = sum(ansMatrix(ansMatrix ~= 0));
            grass_idx = grass_idx + 1;
        end
    elseif contains(file{i}, 'shrub')
        yearStr = extractBetween(file{i}, 'shrub-', '.tif');%提取所需要的年份
        year = str2double(yearStr); 
        if year >= 1990 && year <= 2020
            ansMatrix = readTIF(filepath_t, '灌木', year);
            nonZeroSum(4,shrub_idx) = sum(ansMatrix(ansMatrix ~= 0));
            shrub_idx = shrub_idx + 1;
        end
    elseif contains(file{i}, 'wetland')
        yearStr = extractBetween(file{i}, 'wetland-', '.tif');%提取所需要的年份
        year = str2double(yearStr);
        if year >= 1990 && year <= 2020
            ansMatrix = readTIF(filepath_t, '湿地', year);
            nonZeroSum(5,wetland_idx) = sum(ansMatrix(ansMatrix ~= 0));
            wetland_idx = wetland_idx + 1;
        end
    end
    % switch file{i}
    %     case "cropland*.tif"
    %         nonZeroSum(1,i) = sum(ansMatrix(ansMatrix ~= 0));
    % end  
    fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\覆被图\',file{i},'.jpg');
    %disp(fname);
    if year >= 1990 && year <= 2020
        saveas(gcf, fname);
    end
end

years = 1990 : (1990 + length(nonZeroSum) - 1);
% 创建一个新的图形窗口
figure(1);

% 定义不同的颜色和标记样式
softColors = [
    0.749, 0.682, 0.682;  % 浅灰色
    0.631, 0.722, 0.902;  % 浅蓝色
    0.902, 0.749, 0.722;  % 淡珊瑚色
    0.749, 0.902, 0.749;  % 浅绿色
    0.941, 0.941, 0.627   % 浅黄色
];

set(gcf, 'Color', [0.941, 0.941, 0.941]); % 设置背景为浅灰色

% 绘制每行的数据
for i = 1:size(nonZeroSum, 1)
    % plot(years, nonZeroSum(i, :), ['-', colors{i}], 'LineWidth', 2);
    plot(years, nonZeroSum(i, :), '-', 'Color', softColors(i, :), 'LineWidth', 2);
    hold on; % 保持当前图形以便叠加新的线条
end

% 添加标题和轴标签
title('1990年到2020年国内土地利用率总量的变化');
xlabel('年份');
ylabel('土地利用率总量');
grid on;
set(gca, 'GridAlpha', 0.5, 'GridLineStyle', '--'); % 设置网格线透明度和样式
% 添加图例
legend("田地", "森林" , "草原", "灌木", "湿地", 'Location', 'best');
hold off;
%保存图像
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\覆被图\','线条','.jpg');
disp(fname);
saveas(gcf, fname);

% figure(1);
% plot(years,nonZeroSum(1,:));
% % bar(years,nonZeroSum);
% title("1900年到2020年国内土地利用率总量的变化");
% xlabel("年份");
% ylabel("土地利用率总量");
figure(2);
nonZeroSum_normalized = nonZeroSum./ sum(nonZeroSum(:,1));
b = bar(years,nonZeroSum_normalized',"stacked");

ylim([0 1]);
% 添加标题和标签（可选）
title('土地利用情况(逐年比率)');
xlabel('年份');
ylabel('比率');
land_use_labels = {'田地', '森林', '草原', '灌木', '湿地'}; % 请根据实际情况修改标签
legend(b, land_use_labels, 'Location', 'best'); % 'Location' 参数可以根据需要调整

%保存图像
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\覆被图\','比例','.jpg');
disp(fname);
saveas(gcf, fname);

%% 处理降雨数据
years = 1990 : (1990 + length(nonZeroSum));
%计算1990年和1961年之间的天数差
beginDays = days(datetime('1990-01-01') - datetime('1961-01-01'));
endDays = days(datetime('2020-01-01') - datetime('1961-01-01'));


RainPerYear = zeros(1,(2020 - 1990));
RainPerYearSigma = zeros(1,(2020 - 1990));

d_t = rainData(:,:,1);
NanValueNum = sum(d_t <= 0, 'all');
NanValueNum = 256 * 144 - NanValueNum;%求出有统计值的区块数量

figure(5);
currentY = 1990;
elapsDays = 0;
for i = 1 : (2020 - 1990 + 1)
    stepBeginYear = strcat(num2str(currentY),'-01-01');
    stepEndYear = strcat(num2str(currentY + 1),'-01-01');
    disp(stepBeginYear);
    disp(stepEndYear);

    DAY_step = days(datetime(stepEndYear) - datetime(stepBeginYear));%按照真实年份步进
    disp(DAY_step);  
    elapsDays = elapsDays + DAY_step;

    dataR = zeros(256, 144);
    for j = 1 : DAY_step
        dataR = dataR + rainData(:,:,(beginDays + elapsDays + j));
    end

    dataR(dataR<0) = 0;
    maxDataR = max(max(dataR));
    minDataR = min(min(dataR));
    dataR_normalized = (dataR - minDataR) / (maxDataR - minDataR);

    mesh(dataR');
    colormap("parula"); % 可以选择不同的颜色图
    colorbar;

    titleString = num2str(currentY) + "年降水量";
    currentY = currentY + 1;%年份加1
    title(titleString);
    zlabel('年降水量 单位:(毫米)');

    fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\',stepBeginYear,'.jpg');
    saveas(gcf, fname);

    %年降水量总体
    RainPerYear(i) = sum(dataR, "all");

    %年降水量标准差
    RainPerYearSigma(i) = std(dataR,0,"all");

    pause(0.001);
end

colors = lines(4);
colors(1,:) = [0.5, 0.7, 1];
colors(2,:) = [0.9, 0.6, 0.1];
colors(3,:) = [0.5, 0.8, 0.5];
colors(4,:) = [0.4, 0.6, 0.1];

figure(6);
bar(years, RainPerYear, 'FaceColor', colors(1,:));
title('1990年到2020年国内降水总量的变化情况');
xlabel('年份');
ylabel('降水总量 单位:(毫升)');
grid on;
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\','降水总量','.jpg');
saveas(gcf, fname);

figure(7);
RainPerYear_mean = (RainPerYear./NanValueNum);
plot(years, RainPerYear_mean, 'LineWidth', 2, 'Color', colors(2,:));
title('1990年到2020年国内降水均值的变化情况');
xlabel('年份');
ylabel('年均降水 单位:(毫升)');
grid on;
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\','降水均值','.jpg');
saveas(gcf, fname);


figure(8);
plot(years, RainPerYearSigma, 'LineWidth', 2, 'Color', colors(3,:));
title('1990年到2020年国内降水标准差的变化情况');
xlabel('年份');
ylabel('降水标准差 单位:(毫升)');
grid on;
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\','降水标准差','.jpg');
saveas(gcf, fname);

figure(9);
changeRate_Rain = zeros(1, length(years) - 1);
for i = 2 : length(years)
    changeRate_Rain(i - 1) = (RainPerYear_mean(i) - RainPerYear_mean(i - 1)) / RainPerYear_mean(i - 1);
end
bar(1991:2020, changeRate_Rain, 'LineWidth', 2, 'FaceColor', colors(4,:));
title('1990年到2020年国内降水均值变化率');
xlabel('年份');
ylabel('变化率');
grid on;
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\','降水变化率','.jpg');
saveas(gcf, fname);

%% 问2


