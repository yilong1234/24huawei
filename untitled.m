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

%terrain的图像经纬度范围
terrainInfo = geotiffinfo('C:\Users\Administrator\Desktop\华为杯竞赛\中国数字高程图(1km)\Geo\TIFF\chinadem_geo.tif');
terrainW = terrainInfo.SpatialRef.LongitudeLimits;
terrainN = terrainInfo.SpatialRef.LatitudeLimits;
%气温的经纬度
temperaInfo = geotiffinfo('C:\Users\Administrator\Desktop\华为杯竞赛\日平均数据\1979_avg\19790101_avg.tif');
temperaW = temperaInfo.SpatialRef.LongitudeLimits;
temperaN = temperaInfo.SpatialRef.LatitudeLimits;
% temperaW = [70, 140];
% temperaN = [3, 55];
%降水的经纬度
rainW = [72, 136];
rainN = [18, 54];
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

%第三问需要数据集
croplandData = zeros(78,127,yearsNum);
forestData = zeros(78,127,yearsNum);
grassData = zeros(78,127,yearsNum);
shrubData = zeros(78,127,yearsNum);
wetlandData = zeros(78,127,yearsNum);

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
            type = '农业用地';

            ansMatrix = readTIF(filepath_t);
            

            croplandData(:,:,cropland_idx) = ansMatrix;
            nonZeroSum(1,cropland_idx) = sum(ansMatrix(ansMatrix ~= 0));

            %画出大小144*78的纹路
            totalM = totalM - ansMatrix;

            cropland_idx = cropland_idx + 1;
        end

    elseif contains(file{i}, 'forest')
        yearStr = extractBetween(file{i}, 'forest-', '.tif');%提取所需要的年份
        year = str2double(yearStr);
        if year >= 1990 && year <= 2020
            type = '森林';
            ansMatrix = readTIF(filepath_t);
            forestData(:,:,forest_idx) = ansMatrix;
            nonZeroSum(2,forest_idx) = sum(ansMatrix(ansMatrix ~= 0));
            %画出大小144*78的纹路
            totalM = totalM - ansMatrix;
            
            forest_idx = forest_idx + 1;
        end
    elseif contains(file{i}, 'grass')
        yearStr = extractBetween(file{i}, 'grass-', '.tif');%提取所需要的年份
        year = str2double(yearStr);
        if year >= 1990 && year <= 2020
            type = '草原';
            ansMatrix = readTIF(filepath_t);
            grassData(:,:,grass_idx) = ansMatrix;
            nonZeroSum(3,grass_idx) = sum(ansMatrix(ansMatrix ~= 0));
            %画出大小144*78的纹路
            totalM = totalM - ansMatrix;
            grass_idx = grass_idx + 1;
        end
    elseif contains(file{i}, 'shrub')
        yearStr = extractBetween(file{i}, 'shrub-', '.tif');%提取所需要的年份
        year = str2double(yearStr); 
        if year >= 1990 && year <= 2020
            type = '灌木';
            ansMatrix = readTIF(filepath_t);
            shrubData(:,:,shrub_idx) = ansMatrix;
            nonZeroSum(4,shrub_idx) = sum(ansMatrix(ansMatrix ~= 0));
            %画出大小144*78的纹路
            totalM = totalM - ansMatrix;
            shrub_idx = shrub_idx + 1;
        end
    elseif contains(file{i}, 'wetland')
        yearStr = extractBetween(file{i}, 'wetland-', '.tif');%提取所需要的年份
        year = str2double(yearStr);
        if year >= 1990 && year <= 2020
            type = '湿地';
            ansMatrix = readTIF(filepath_t);
            wetlandData(:,:,wetland_idx) = ansMatrix;
            nonZeroSum(5,wetland_idx) = sum(ansMatrix(ansMatrix ~= 0));
            %画出大小144*78的纹路
            totalM = totalM - ansMatrix;
            wetland_idx = wetland_idx + 1;
        end
    end

    if(year >= 1990 && year <= 2020)
        pcolor(ansMatrix);
        % 设置坐标轴刻度
        axis ij;
        axis tight;
        % 移除网格线
        clim([0 1]); % 设置颜色条的最小值和最大值
        set(gca, 'GridLineStyle', 'none'); % 关闭网格线

        % 设置标题和其他属性
        titleValue = strcat('土地利用情况：',type);
        title(titleValue);
        xlabelStr = strcat('当前年份为', num2str(year));
        xlabel(xlabelStr);

        % 显示颜色条
        colormap(jet); % 可以选择不同的颜色图
        colorbar;

        % fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\覆被图\',file{i},'.jpg');
        % saveas(gcf, fname);
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
% fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\覆被图\','线条','.jpg');
% disp(fname);
% saveas(gcf, fname);

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
% fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\覆被图\','比例','.jpg');
% disp(fname);
% saveas(gcf, fname);

%% 处理降雨数据
years = 1990 : (1990 + length(nonZeroSum));
%计算1990年和1961年之间的天数差
beginDays = days(datetime('1990-01-01') - datetime('1961-01-01'));
endDays = days(datetime('2020-01-01') - datetime('1961-01-01'));

RainPerYear = zeros(1,(2020 - 1990));
RainPerYearSigma = zeros(1,(2020 - 1990));
torrentialrain = zeros(256, 144,(2020 - 1990 + 1));

dataR = zeros(256, 144, (2020 - 1990 + 1));

d_t = rainData(:,:,1);
NanValueNum = sum(d_t <= 0, 'all');
NanValueNum = 256 * 144 - NanValueNum;%求出有统计值的区块数量

%画出轮廓
d_t(d_t > 0) = 0;


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

    
    for j = 1 : DAY_step
         %当前雨天
        currentRainData = rainData(:, :, (beginDays + elapsDays + j));
        currentRainData(currentRainData < 0) = 0;

        dataR(:,:,i) = dataR(:,:,i) + currentRainData;

        %找出索引
        idx = find(currentRainData > 50);
        [row, col] = ind2sub(size(currentRainData), idx);
        %记录暴雨数
        torrentialrain(row, col, i) = torrentialrain(row, col, i) + 1;
    end
    
    
    % maxDataR = max(max(dataR));
    % minDataR = min(min(dataR));
    % dataR_normalized = (dataR - minDataR) / (maxDataR - minDataR);
    
    % mesh(dataR(:,:,i)');

    rainstepN = (rainN(2) - rainN(1))/144;
    rainstepW = (rainW(2) - rainW(1))/256;
    x = (rainN(1) : rainstepN : rainN(2) - rainstepN);
    y = (rainW(1) : rainstepW : rainW(2) - rainstepW);

    currentRain = dataR(:,:,i)';
    currentRain(currentRain == 0) = 0;
    pause(1);
    mesh(y, x, currentRain);

    colormap("parula"); % 可以选择不同的颜色图
    colorbar;
    view(0,90);
    zlabel('年降水量 单位:(毫米)');
    xlabel('经度');
    ylabel('维度');
    titleString = num2str(currentY) + "年降水量";
    currentY = currentY + 1;%年份加1
    title(titleString);

    % fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\',stepBeginYear,'.jpg');
    % saveas(gcf, fname);

    %年降水量总体
    RainPerYear(i) = sum(dataR(:,:,i), "all");

    %年降水量标准差
    RainPerYearSigma(i) = std(dataR(:,:,i),0,"all");
end

rainMean = zeros(256, 144);
for i = 1 : (2020 - 1990 + 1)
    rainMean = rainMean + dataR(:,:,i);
end
mesh(fliplr(rainMean./(2020 - 1990 + 1)));
title("不同地区年均降水量")
ylabel('降水总量 单位:(毫升)');
colorbar;
% fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\','不同地区年均降水量','.jpg');
% saveas(gcf, fname);

%求出暴雨率,暴雨天数/365
% torrentialrain = torrentialrain./365;
mesh(torrentialrain(:, :, 1)');

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

% fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\','降水总量','.jpg');
% saveas(gcf, fname);

figure(7);
RainPerYear_mean = (RainPerYear./NanValueNum);
plot(years, RainPerYear_mean, 'LineWidth', 2, 'Color', colors(2,:));
title('1990年到2020年国内降水均值的变化情况');
xlabel('年份');
ylabel('年均降水 单位:(毫升)');
grid on;
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
% fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\','降水均值','.jpg');
% saveas(gcf, fname);


figure(8);
plot(years, RainPerYearSigma, 'LineWidth', 2, 'Color', colors(3,:));
title('1990年到2020年国内降水标准差的变化情况');
xlabel('年份');
ylabel('降水标准差 单位:(毫升)');
grid on;
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
% fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\','降水标准差','.jpg');
% saveas(gcf, fname);

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

figure(13);
totalChangeRate_Rain = zeros(1, length(years) - 1);
for i = 2 : length(years)
    totalChangeRate_Rain(i - 1) = (RainPerYear(i) - RainPerYear(i - 1)) / RainPerYear(i - 1);
end
bar(1991:2020, totalChangeRate_Rain, 'LineWidth', 2, 'FaceColor', colors(4,:));
title('1990年到2020年国内降水总量变化率');
xlabel('年份');
ylabel('变化率');
grid on;
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
% fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问1\降水量\','降水变化率','.jpg');
% saveas(gcf, fname);

%% 问2 读取地形数据
terrain = readTIF('C:\Users\Administrator\Desktop\华为杯竞赛\中国数字高程图(1km)\Geo\TIFF\chinadem_geo.tif');

% 目标尺寸
targetHeight = 256;
targetWidth = 144;
% 使用最大池化改变数组尺寸，可以单独作为理论小节
resizedTerrain = MaxPool(terrain, targetHeight, targetWidth);
resizedTerrain(resizedTerrain == -32768) = 0;

figure(10);
x_terrain_stepN = (terrainN(2) - terrainN(1))/256;
x_terrain_stepW = (terrainW(2) - terrainW(1))/144;

x = (terrainW(1) : x_terrain_stepW : terrainW(2) - x_terrain_stepW);
y = (terrainN(2) : -x_terrain_stepN : terrainN(1) + x_terrain_stepN);

mesh(x,y,resizedTerrain);
title('最大池化后的地形结果');
zlabel('高度 单位(米)');
xlabel('纬度 单位(度)');
ylabel('经度 单位(度)');
colorbar;
view(0,70);
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问2\地形数据\','池化后的地形图','.jpg');
saveas(gcf, fname);

%% 读取气温数据 1990 -2018
figure(11);

temperation_mean = zeros(400, 700,(2018 - 1990 + 1));
temperation_mean_resized = zeros(targetHeight, targetWidth,(2018 - 1990 + 1));
currentY = 1990;


for i = 1 : (2018 - 1990 + 1)
    pause(0.001);

    temper_str = dir(strcat('C:\Users\Administrator\Desktop\华为杯竞赛\日平均数据\', num2str(1990 + i - 1),'_avg'));
    temper_folder = string({temper_str.folder}');
    temper_file = string({temper_str.name}');
    temper_filepath = strcat(temper_folder, '\', temper_file);

    temperation = zeros(400, 700);

    for j = 3 : length(temper_file)
        filepath_t = temper_filepath{j};
        if(contains(filepath_t, 'tif'))
            % disp(filepath_t);
            temper_Matrix = readTIF(filepath_t);
            temperation = temperation + temper_Matrix;%求和           
        end     
    end
    stepBeginYear = strcat(num2str(currentY),'-01-01');
    stepEndYear = strcat(num2str(currentY + 1),'-01-01');
    DAY_step = days(datetime(stepEndYear) - datetime(stepBeginYear));%按照真实年份步进

    temperation_mean(:,:,i) = temperation./DAY_step;
    temperation_mean(temperation_mean <= -111) = -111;
    % temperation_mean_resize = MaxPool(temperation_mean, targetHeight, targetWidth);
    % temperation_mean_resized(:,:,i) = imresize(temperation_mean(:,:,i), [targetHeight, targetWidth], 'bicubic');

    x_temper_stepN = (temperaN(2) - temperaN(1))/400;
    x_temper_stepW = (temperaW(2) - temperaW(1))/700;

    y = (temperaW(1) : x_temper_stepW : temperaW(2) - x_temper_stepW);
    x = (temperaN(2) : -x_temper_stepN : temperaN(1) + x_temper_stepN);

    surf(y,x,temperation_mean(:,:,i), 'EdgeColor', 'none');
    % mesh(temperation_mean_resize');
    title(strcat(num2str(1990 + i - 1), '年平均温度'));
    zlabel('温度 单位(摄氏度)');
    xlabel('经度方向');
    ylabel('纬度方向');
    view(0, 90);
    colorbar;
    colormap(jet);

    fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问2\地形数据\',strcat(num2str(1990 + i - 1), '年平均温度'),'.jpg');
    saveas(gcf, fname);

    currentY = currentY + 1;
end
%% 多元线性拟合
%以海拔、降水、气温构建回归方程
%取一年的数据构建方程
rain1990 = dataR(:,:,1);
rain1990 = fliplr(rain1990);
tor_rain = torrentialrain(:,:,1) + d_t; 
tor_rain = fliplr(tor_rain);
temper1990 = temperation_mean(:,:,1)';


mesh(tor_rain);
mesh(rain1990);
mesh(temper1990);

terrain = terrain';
mesh(terrain);




% xySize = 2000;
% XandY = zeros(6, xySize);

terrainW_step = (terrainW(2) - terrainW(1))/10892;
terrainN_step = (terrainN(2) - terrainN(1))/5385;

% 线性索引
idx = find(tor_rain>=0);
[row, col] = ind2sub(size(tor_rain), idx);
% row为经度方向，col为纬度方向
% 计算所处位置的经纬度
XandY = zeros(6, length(row));
%在纬度相同的情况下每隔1度，距离相差约100千米；
stepW = 100;
%在经度度相同的情况下每隔1度，距离相差约111.32千米；
stepN = 111.32;


for i = 1 : length(row)
    XandY(1, i) = row(i);
    XandY(2, i) = col(i);

    % 计算所处位置的经纬度
    currentW = rainW(2) - ((256 - (row(i))) * 0.25);
    currentN = rainN(2) - (col(i) * 0.25);


    %计算在温度数据集下的坐标
    % temperaRow = temperaW(2) - currentW;
    % temperaRow = temperaRow/0.1;
    temperaRow = (currentW - temperaW(1)) / 0.1;
    temperaCol = temperaN(2) - currentN;
    temperaCol = temperaCol/0.1;


    % XandY(3, i) = temper1990(round(temperaRow), round(temperaCol));
    % 确保中心点在有效范围内

    if temperaRow >= 3 && temperaRow <= (size(temper1990, 1) - 2) && temperaCol >= 3 && temperaCol <= (size(temper1990, 2) - 2)
        % 提取 5x5 子矩阵
        subTemper1990 = temper1990(floor(temperaRow) - 2:floor(temperaRow) + 2, floor(temperaCol) - 2:floor(temperaCol) + 2);
        nonZeroValues = subTemper1990(abs(subTemper1990) > 0);
        if size(nonZeroValues)>0
            % 将非零值存入 XandY 的第三个元素作为细胞数组
            XandY(3, i) = max(nonZeroValues(:));
        else
            XandY(3, i) = 0;
        end

    end


    %计算在地形数据集下的坐标
    terrainRow = (currentW - terrainW(1)) / terrainW_step;
    terrainCol = (terrainN(2) - currentN) / terrainN_step;
    if terrainRow >= 3 && terrainRow <= (size(terrain, 1) - 2) && terrainCol >= 3 && terrainCol <= (size(terrain, 2) - 2)
        % 提取 9x9 子矩阵
        subTerrain = terrain(floor(terrainRow) - 2:floor(terrainRow) + 2, floor(terrainCol) - 2:floor(terrainCol) + 2);

        % 找到子矩阵中的最大值
        XandY(4, i) = max(subTerrain(:));
    end
    % 计算固定纬度（即经度方向）的距离
    % XandY(4, i) = terrain(floor(terrainRow), floor(terrainCol));
    %存放暴雨数和降雨量
    XandY(5, i) = rain1990(row(i), col(i));
    XandY(6, i) = tor_rain(row(i), col(i));
    
    if i == 1
        disp("i = " + i);

        disp(currentW);
        disp(currentN); 
        disp(row(i));
        disp(col(i));

        disp(temperaRow);
        disp(temperaCol);

        disp(terrainW(2));
        disp(terrainN(2));

        disp(terrainRow);
        disp(terrainCol);
    end
    
end

%% 去除错误的采样数据部分开始拟合
%错误的采样
threshold = -32768;
errIdx = zeros(1,length(row));
idx_t = 1;
for i = 1 : length(row)
    if XandY(4, i) == threshold
        errIdx(idx_t) = i;
        idx_t = idx_t + 1;
    end
end
errIdx(errIdx==0) = [];
for i = length(errIdx) :-1: 1 
    XandY(:,[errIdx(i)]) = [];
end
%拟合
preRegress = XandY(3:6, :);
preRegress = preRegress';
X = all2one(preRegress(:,1:3)); %前6列为自变量
Y = preRegress(:,4);   %最后一列为因变量
result = fitlm(X, Y);
disp(result);
%% 建立相似矩阵分析各个参数之间的关系
correlationMatrix = corrcoef(X);
% disp(correlationMatrix);
figure;
imagesc(correlationMatrix);
colorbar; % 添加颜色条
title('特征间的相关矩阵');
% xlabel('特征编号');
% ylabel('特征编号');

% 设置 x 轴和 y 轴标签为特征名称（如果有的话）
set(gca, 'XTick', 1:3, 'YTick', 1:3, ...
     'XTickLabel', {'气温', '海拔', '降雨'}, ...
     'YTickLabel', {'气温', '海拔', '降雨'});
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问2\地形数据\','气温,海拔降雨的相关矩阵','.jpg');
saveas(gcf, fname);

%% 绘制模型验证相关
plotResiduals(result, 'fitted'); % 绘制拟合值与残差的关系图
plotResiduals(result, 'probability'); % 绘制残差的正态概率图
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问2\地形数据\','残差的正态概率图','.jpg');
saveas(gcf, fname);
plotResiduals(result, 'fitted'); % 观察残差随拟合值的变化情况
%% 使用2018年数据验证
rain2018 = dataR(:,:,end);
rain2018 = fliplr(rain2018);
temper2018 = temperation_mean(:,:,(2018 - 1990))';
tor_rain2018  = torrentialrain(:,:,(2018 - 1990)) + d_t; 
tor_rain2018  = fliplr(tor_rain2018);

random_idx = randperm(length(idx), 100);

XandY2018 = zeros(6, length(row));
for i = 1 : length(row)
    XandY2018(1, i) = row(i);
    XandY2018(2, i) = col(i);

    % 计算所处位置的经纬度
    currentW = rainW(2) - ((256 - (row(i))) * 0.25);
    currentN = rainN(2) - (col(i) * 0.25);

    %计算在温度数据集下的坐标
    % temperaRow = temperaW(2) - currentW;
    % temperaRow = temperaRow/0.1;
    temperaRow = (currentW - temperaW(1)) / 0.1;
    temperaCol = temperaN(2) - currentN;
    temperaCol = temperaCol/0.1;

    % XandY(3, i) = temper1990(round(temperaRow), round(temperaCol));
    % 确保中心点在有效范围内

    if temperaRow >= 3 && temperaRow <= (size(temper2018, 1) - 2) && temperaCol >= 3 && temperaCol <= (size(temper2018, 2) - 2)
        % 提取 5x5 子矩阵
        subTemper1990 = temper2018(floor(temperaRow) - 2:floor(temperaRow) + 2, floor(temperaCol) - 2:floor(temperaCol) + 2);
        nonZeroValues = subTemper1990(abs(subTemper1990) > 0);
        if size(nonZeroValues)>0
            % 将非零值存入 XandY 的第三个元素作为细胞数组
            XandY2018(3, i) = max(nonZeroValues(:));
        else
            XandY2018(3, i) = 0;
        end

    end
    
    %计算在地形数据集下的坐标
    terrainRow = (currentW - terrainW(1)) / terrainW_step;
    terrainCol = (terrainN(2) - currentN) / terrainN_step;
    if terrainRow >= 3 && terrainRow <= (size(terrain, 1) - 2) && terrainCol >= 3 && terrainCol <= (size(terrain, 2) - 2)
        % 提取 9x9 子矩阵
        subTerrain = terrain(floor(terrainRow) - 2:floor(terrainRow) + 2, floor(terrainCol) - 2:floor(terrainCol) + 2);

        % 找到子矩阵中的最大值
        XandY2018(4, i) = max(subTerrain(:));
    end
    % 计算固定纬度（即经度方向）的距离
    % XandY(4, i) = terrain(floor(terrainRow), floor(terrainCol));
    %存放暴雨数和降雨量
    XandY2018(5, i) = rain2018(row(i), col(i));
    XandY2018(6, i) = tor_rain2018(row(i), col(i));
    
end
selected_elements = XandY2018(:,random_idx);
threshold = -32768;
errIdx = zeros(1,length(100));
idx_t = 1;
for i = 1 : 100
    if selected_elements(4, i) == threshold
        errIdx(idx_t) = i;
        idx_t = idx_t + 1;
    end
end
errIdx(errIdx==0) = [];
if ~isempty(errIdx)
    for i = length(errIdx)
        selected_elements(:,[errIdx(i)]) = [];
    end
end
selected_elements = selected_elements';
X = selected_elements(:,3:5); %前6列为自变量
Y = selected_elements(:,6);
y_predict = predict(result, X);
% for i = 1 : length(Y)
%     y_predict(i) = 0.238 * X(i,1) + 1.8918e-5 * X(i,2) + 0.0044929 * X(i,3);
% end
Y = Y';
figure(12);
set(gcf, 'WindowState', 'maximized');
plot(Y,'o-');
hold on;
plot(y_predict,'*-');
legend('真实值','预测值');
hold off;
title('模型对2018年极端暴雨情况的验证');
ylabel('年暴雨总次数');
xlabel('不同采样点');
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问2\地形数据\','拟合的结果','.jpg');
saveas(gcf, fname);


%% 问3 提取训练集
%建立不同经纬度下的映射关系，以土地利用率为主要数组
%土地利用率的经纬度,主要的经纬度
usedLandInfo = geotiffinfo('C:\Users\Administrator\Desktop\华为杯竞赛\中国大陆0.5°土地利用和覆盖变化数据集(1900-2019年)\数据实体\cropland-1900.tif');
usedLandW = usedLandInfo.SpatialRef.LongitudeLimits;
usedLandN = usedLandInfo.SpatialRef.LatitudeLimits;
%terrain的图像经纬度范围
terrainInfo = geotiffinfo('C:\Users\Administrator\Desktop\华为杯竞赛\中国数字高程图(1km)\Geo\TIFF\chinadem_geo.tif');
terrainW = terrainInfo.SpatialRef.LongitudeLimits;
terrainN = terrainInfo.SpatialRef.LatitudeLimits;
%气温的经纬度
temperaInfo = geotiffinfo('C:\Users\Administrator\Desktop\华为杯竞赛\日平均数据\1979_avg\19790101_avg.tif');
temperaW = temperaInfo.SpatialRef.LongitudeLimits;
temperaN = temperaInfo.SpatialRef.LatitudeLimits;
%降水的经纬度
rainW = [72, 136];
rainN = [18, 54];


%初始化训练集比例为2/3
years = 29;
SVMpredictionX = zeros(years * 4184,8);
SVMY = zeros(years * 4184,1);

evl_num_total = croplandData(:,:,1) + forestData(:,:,1) + grassData(:,:,1) + shrubData(:,:,1) + wetlandData(:,:,1);
%转置统一数据
evl_num_total_T = evl_num_total';
idx_usedLand = find(evl_num_total_T > 0);
[row_usedLand, col_usedLand] = ind2sub(size(evl_num_total_T), idx_usedLand);


% mesh(evl_num_total_T);
% mesh(temper2018);
for i = 1 : years
    for j = 1 : 4184

        % 计算所处位置的经纬度
        currentW = usedLandW(2) - ((127 - (row_usedLand(j))) * 0.5);
        currentN = usedLandN(2) - (col_usedLand(j) * 0.5);
        % 
        % 
        % %拿到温度
        temperaRow = round((currentW - temperaW(1)) / 0.1);
        temperaCol = round((temperaN(2) - currentN) / 0.1);
        temperationT = temperation_mean(:,:,i)';
        SVMpredictionX((i - 1) * 4184 + j, 1) = temperationT(temperaRow, temperaCol);%第1个参数代表温度

        % %拿到降水
        raincol = round((currentW - rainW(1)) / 0.25);
        rainRow = round((currentN - rainN(1)) / 0.25);
        rainT = dataR(:,:,i)';
        % mesh(rainT);
        if(rainRow == 0)
            rainRow = 1;
        end
        SVMpredictionX((i - 1) * 4184 + j, 3) = rainT(rainRow, raincol);%第3个参数代表降水

        %拿到海拔
        terrainRow = round((currentW - terrainW(1)) / terrainW_step);
        terrainCol = round((terrainN(2) - currentN) / terrainN_step);
        SVMpredictionX((i - 1) * 4184 + j, 2) = terrain(terrainRow, terrainCol);%第2个参数代表海拔

        %拿到土地利用率
        croplandDataT = croplandData(:,:, i)';
        forestDataT = forestData(:,:, i)';
        grassDataT = grassData(:,:, i)';
        shrubDataT = shrubData(:,:, i)';
        wetlandDataT = wetlandData(:,:, i)';

        SVMpredictionX((i - 1) * 4184 + j, 4) = croplandDataT(row_usedLand(j), col_usedLand(j));%第4个参数代表农业用地
        SVMpredictionX((i - 1) * 4184 + j, 5) = forestDataT(row_usedLand(j), col_usedLand(j));%第5个参数代表森林
        SVMpredictionX((i - 1) * 4184 + j, 6) = grassDataT(row_usedLand(j), col_usedLand(j));%第6个参数代表草原
        SVMpredictionX((i - 1) * 4184 + j, 7) = shrubDataT(row_usedLand(j), col_usedLand(j));%第7个参数代表灌木
        SVMpredictionX((i - 1) * 4184 + j, 8) = wetlandDataT(row_usedLand(j), col_usedLand(j));%第8个参数代表湿地

        %拿到暴雨数，并且以年均暴雨5次及以上为1，不超过5次的为0
        torrentialrainT = torrentialrain(:,:,i)';
        SVMY((i - 1) * 4184 + j, 1) = torrentialrainT(rainRow, raincol);
        % mesh(torrentialrainT + d_t');
        % if(torrentialrainT(rainRow, raincol) >= 10)
        %     SVMY((i - 1) * 4184 + j, 1) = 1;
        % else
        %     SVMY((i - 1) * 4184 + j, 1) = 0;
        % end

        if j == 1 && i == 1
            disp(i);
            disp(row_usedLand(j));
            disp(col_usedLand(j));
            disp(currentW);
            disp(currentN);

            disp(temperaRow);
            disp(temperaCol);

            disp(rainRow);
            disp(raincol);

            disp(terrainRow);
            disp(terrainCol);
        end
        

    end
end
% %对异常值做剔除，用线性插值方法
% col_needProcess = 1:3; % 处理所有列
% SVMX = SVMpredictionX;
% % 对每一列分别进行异常值处理
% for col = col_needProcess
%     % 使用 filloutliers 函数按列方向处理异常值
%     SVMX(:, col) = filloutliers(SVMX(:, col),'linear');
% end
%手动清洗数据
%错误数据的清理
SVMX = SVMpredictionX;
errIdx = zeros(1,years * 4184);
idx_t = 1;
for i = 1 : (years * 4184)
    if SVMX(i, 1) == -111 || SVMX(i, 2) == -32768
        errIdx(idx_t) = i;
        idx_t = idx_t + 1;
    end
end
errIdx(errIdx==0) = [];
for i = length(errIdx) :-1: 1
    SVMX([errIdx(i)],:) = [];
    SVMY([errIdx(i)],:) = [];
end

div = 4184 - length(errIdx) / years;
trainNum = div * 0.7;
TestNum = div * 0.3;
%训练集
SVMX_train = SVMX(1: div, 2 : 8);
SVMY_train = SVMY(1: div, :);
%测试集
SVMX_Test = SVMX(div + trainNum + 1 : 2 * div, 2 : 8);
SVMY_Test = SVMY(div + trainNum + 1 : 2 * div, :);
%训练集归一化
SVMX_train(:, 1:2) = all2one(SVMX_train(:, 1:2));
SVMX_Test(:, 1:2) = all2one(SVMX_Test(:, 1:2));
%% 相关矩阵
correlationMatrix = corrcoef(SVMX_train);
% disp(correlationMatrix);
figure;
imagesc(correlationMatrix);
colorbar; % 添加颜色条
title('特征间的相关矩阵');
% xlabel('特征编号');
% ylabel('特征编号');

% 设置 x 轴和 y 轴标签为特征名称（如果有的话）
set(gca, 'XTick', 1:8, 'YTick', 1:8, ...
     'XTickLabel', {'气温', '海拔', '降雨','农业用地','森林','草地','灌木','湿地'}, ...
     'YTickLabel', {'气温', '海拔', '降雨','农业用地','森林','草地','灌木','湿地'});
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问3\','问3模型参数的相关矩阵','.jpg');
saveas(gcf, fname);
%% 训练模型

svmModel = fitrsvm(SVMX_train, SVMY_train,'KernelFunction','linear');
disp(svmModel);
predictedY = predict(svmModel, SVMX_Test);
%% 验证结果
figure;
residuals = SVMY_Test - predictedY;
scatter(predictedY, residuals, 'filled');
xlabel('预测值');
ylabel('残差');
title('模型拟合的残差图');
grid on;
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问3\','模型拟合的残差图','.jpg');
saveas(gcf, fname);

figure;
set(gcf, 'WindowState', 'maximized');
scatter(SVMY_Test, predictedY, 'filled');
hold on;
plot(min(SVMY_Test):max(SVMY_Test), min(SVMY_Test):max(SVMY_Test), 'r--');
xlabel('实际值');
ylabel('预测值');
title('实际值 vs. 预测值');
legend('数据点', '理想线');
grid on;
hold off;
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问3\','实际值 vs. 预测值','.jpg');
saveas(gcf, fname);

% 计算模型的准确性
mse = mean((predictedY - SVMY_Test).^2); % 均方误差
rmse = sqrt(mse); % 均方根误差
rsquared = 1 - sum((SVMY_Test - predictedY).^2) / sum((SVMY_Test - mean(SVMY_Test)).^2); % 决定系数

% 显示结果
fprintf('均方误差 (MSE): %.4f\n', mse);
fprintf('均方根误差 (RMSE): %.4f\n', rmse);
fprintf('决定系数 (R²): %.4f\n', rsquared);

figure(19);
set(gcf, 'WindowState', 'maximized');
plot(predictedY,'-o');
hold on;
plot(SVMY_Test,'-*');
legend('预测值','真实值');
hold off;
title('模型对测试集的验证');
ylabel('年暴雨总次数');
xlabel('不同采样点');
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问3\','拟合的结果','.jpg');
saveas(gcf, fname);



%% 准备2025年到2035年的数据
% %求出每年不同覆被的变化率，以便于预测2025到2035的数据
% mesh(croplandData(:,:,30));
% nonZeroSum_predict = zeros(5, 30);
% for i = 1 : 5
%     for j = 2 : 30
%         nonZeroSum_predict(i, j - 1) = (nonZeroSum_normalized(i, j) - nonZeroSum_normalized(i, j - 1)) / nonZeroSum_normalized(i, j - 1);
%     end
% end
% nonZeroSum_predict = nonZeroSum_predict(:,1 : 29);
% plot(nonZeroSum_predict(1, :));
% hold on;
% plot(nonZeroSum_predict(2, :));
% plot(nonZeroSum_predict(3, :));
% plot(nonZeroSum_predict(4, :));
% plot(nonZeroSum_predict(5, :));
% hold off;
% legend("田地", "森林" , "草原", "灌木", "湿地", 'Location', 'best');

%计算2025-2035的土地利用情况
% 定义原始的和新的索引
original_indices = 1:30;  % 原始的30个点的索引
new_indices = linspace(1, 30, 45);  % 新的45个点的索引，包括原始点和插值点

% 初始化一个新的空矩阵来存储结果
cropData_B = zeros(78, 127, 45);
forestData_B = zeros(78, 127, 45);
grassData_B = zeros(78, 127, 45);
shrubData_B = zeros(78, 127, 45);
wetlandData_B = zeros(78, 127, 45);

%% 预测土地利用率等数据
% [X, Y, T] = size(croplandData);
% numForecastSteps = 15;
% croplandData_forecasts = zeros(X, Y, numForecastSteps); % 假设 numForecastSteps 是你想要预测的时间步数
% forestData_forecasts = zeros(X, Y, numForecastSteps); % 假设 numForecastSteps 是你想要预测的时间步数
% grassData_forecasts = zeros(X, Y, numForecastSteps); % 假设 numForecastSteps 是你想要预测的时间步数
% shrubData_forecasts = zeros(X, Y, numForecastSteps); % 假设 numForecastSteps 是你想要预测的时间步数
% wetlandData_forecasts = zeros(X, Y, numForecastSteps); % 假设 numForecastSteps 是你想要预测的时间步数
% 
% for x = 1:X
%     for y = 1:Y
%         timeSeries = croplandData(x, y, :); % 提取当前 (x, y) 位置的时间序列
%         timeSeries = timeSeries(:);
%         if ~all(isnan(timeSeries)) % 确保时间序列不全是 NaN
%             EstMdl = estimate(Mdl, timeSeries); % 拟合 ARIMA 模型
%             forecastedValues = forecast(EstMdl, numForecastSteps, 'Y0', timeSeries); % 进行预测
%             croplandData_forecasts(x, y, :) = forecastedValues; % 存储预测结果
%         end
%     end
% end
% disp(EstMdl);
% 
% for x = 1:X
%     for y = 1:Y
%         timeSeries = forestData(x, y, :); % 提取当前 (x, y) 位置的时间序列
%         timeSeries = timeSeries(:);
%         if ~all(isnan(timeSeries)) % 确保时间序列不全是 NaN
%             EstMdl = estimate(Mdl, timeSeries); % 拟合 ARIMA 模型
%             forecastedValues = forecast(EstMdl, numForecastSteps, 'Y0', timeSeries); % 进行预测
%             forestData_forecasts(x, y, :) = forecastedValues; % 存储预测结果
%         end
%     end
% end
% 
% for x = 1:X
%     for y = 1:Y
%         timeSeries = grassData(x, y, :); % 提取当前 (x, y) 位置的时间序列
%         timeSeries = timeSeries(:);
%         if ~all(isnan(timeSeries)) % 确保时间序列不全是 NaN
%             EstMdl = estimate(Mdl, timeSeries); % 拟合 ARIMA 模型
%             forecastedValues = forecast(EstMdl, numForecastSteps, 'Y0', timeSeries); % 进行预测
%             grassData_forecasts(x, y, :) = forecastedValues; % 存储预测结果
%         end
%     end
% end
% 
% for x = 1:X
%     for y = 1:Y
%         timeSeries = shrubData(x, y, :); % 提取当前 (x, y) 位置的时间序列
%         timeSeries = timeSeries(:);
%         if ~all(isnan(timeSeries)) % 确保时间序列不全是 NaN
%             EstMdl = estimate(Mdl, timeSeries); % 拟合 ARIMA 模型
%             forecastedValues = forecast(EstMdl, numForecastSteps, 'Y0', timeSeries); % 进行预测
%             shrubData_forecasts(x, y, :) = forecastedValues; % 存储预测结果
%         end
%     end
% end
% 
% for x = 1:X
%     for y = 1:Y
%         timeSeries = wetlandData(x, y, :); % 提取当前 (x, y) 位置的时间序列
%         timeSeries = timeSeries(:);
%         if ~all(isnan(timeSeries)) % 确保时间序列不全是 NaN
%             EstMdl = estimate(Mdl, timeSeries); % 拟合 ARIMA 模型
%             forecastedValues = forecast(EstMdl, numForecastSteps, 'Y0', timeSeries); % 进行预测
%             wetlandData_forecasts(x, y, :) = forecastedValues; % 存储预测结果
%         end
%     end
% end


%%
for i = 1:78
    for j = 1:127
        cropData_B(i, j, :) = interp1(original_indices, squeeze(croplandData(i, j, :)), new_indices, 'linear');
        forestData_B(i, j, :) = interp1(original_indices, squeeze(forestData(i, j, :)), new_indices, 'linear');
        grassData_B(i, j, :) = interp1(original_indices, squeeze(grassData(i, j, :)), new_indices, 'linear');
        shrubData_B(i, j, :) = interp1(original_indices, squeeze(shrubData(i, j, :)), new_indices, 'linear');
        wetlandData_B(i, j, :) = interp1(original_indices, squeeze(wetlandData(i, j, :)), new_indices, 'linear');
    end
end

cropData_B = cropData_B(:,:, 35:45);
forestData_B = forestData_B(:,:, 35:45);
grassData_B = grassData_B(:,:, 35:45);
shrubData_B = shrubData_B(:,:, 35:45);
wetlandData_B = wetlandData_B(:,:, 35:45);

cropData_C = zeros(127,78);
forestData_C = zeros(127,78);
grassData_C = zeros(127,78);
shrubData_C = zeros(127,78);
wetlandData_C = zeros(127,78);
for i = 1 : 11
    cropData_C(:,:,i) = cropData_B(:,:,i)';
    forestData_C(:,:,i) = forestData_B(:,:,i)';
    grassData_C(:,:,i) = grassData_B(:,:,i)';
    shrubData_C(:,:,i) = shrubData_B(:,:,i)';
    wetlandData_C(:,:,i) = wetlandData_B(:,:,i)';
end
mesh(cropData_C(:,:,1));

dataTotal = zeros(256, 144);
%降水具有时空变异性和不可控性，故取30年均值
for i = 1 : (2020 - 1990 + 1)
    dataTotal = dataTotal + dataR(:,:,i);
end
dataTotal = dataTotal./(2020 - 1990 + 1);
dataTotal = dataTotal';
mesh(dataTotal);
%%
%2025-2035不同地区的暴雨数量
torRain25_35 = zeros(127,78,11);

preparetoPerdict_partA = zeros(4184,2);

%提取预测数据集
for i = 1 : 4184
    currentW = usedLandW(2) - ((127 - (row_usedLand(i))) * 0.5);
    currentN = usedLandN(2) - (col_usedLand(i) * 0.5);

    %拿到降水
    raincol = round((currentW - rainW(1)) / 0.25);
    rainRow = round((currentN - rainN(1)) / 0.25);

    % mesh(dataTotal);
    if(rainRow == 0)
        rainRow = 1;
    end
    preparetoPerdict_partA(i,2) = dataTotal(rainRow, raincol);
    

    %拿到海拔，取归一化值
    terrainRow = round((currentW - terrainW(1)) / terrainW_step);
    terrainCol = round((terrainN(2) - currentN) / terrainN_step);
    preparetoPerdict_partA(i,1) = terrain(terrainRow, terrainCol);%第2个参数代表海拔

end


usedLandDATA = zeros(4184, 7, 11);
for i = 1 : 11
    for j = 1 : 4184

        %拿到土地利用率

        usedLandDATA(j, 1, i) = cropData_C(row_usedLand(j),col_usedLand(j), i);
        usedLandDATA(j, 2, i) = forestData_C(row_usedLand(j),col_usedLand(j), i);
        usedLandDATA(j, 3, i) = grassData_C(row_usedLand(j),col_usedLand(j), i);
        usedLandDATA(j, 4, i) = shrubData_C(row_usedLand(j),col_usedLand(j), i);
        usedLandDATA(j, 5, i) = wetlandData_C(row_usedLand(j),col_usedLand(j), i);
        usedLandDATA(j, 6, i) = row_usedLand(j);
        usedLandDATA(j, 7, i) = col_usedLand(j);
        if j == 1 && i == 1
            disp(i);
            disp(row_usedLand(j));
            disp(col_usedLand(j));
        end
    end
end


for i = 1 : 11

    X_p = [preparetoPerdict_partA, usedLandDATA(:,:,i)];

    %数据清洗
    errIdx = zeros(1,4184);
    idx_t = 1;
    for j = 1 : 4184
        if X_p(j, 1) == -32768
            errIdx(idx_t) = j;
            idx_t = idx_t + 1;
        end
    end
    errIdx(errIdx==0) = [];
    for j = length(errIdx) :-1: 1
        X_p([errIdx(j)],:) = [];
        X_p([errIdx(j)],:) = [];
    end

    row_col = X_p(:,8 : 9);
    X_p = X_p(:,1 : 7);

    X_p(:, 1:2) = all2one(X_p(:, 1:2));
    X_p(:, 1:2) = all2one(X_p(:, 1:2));



    predictedY = predict(svmModel, X_p);

    rain_mat = zeros(127,78);
    for m = 1 : 3890
        rain_mat(row_col(m, 1), row_col(m, 2)) =  predictedY(m);
    end
    
    pause(1);
    x = 55.5751 : -0.5 : 16.5750; % 北纬从高到低
    y = 72.2250 : 0.5 : 135.7250;  % 经度从西向东
    [X, Y] = meshgrid(y, x);
    X = X';
    Y = Y';


    mesh(X(1:127, 1:78), Y(1:127, 1:78), rain_mat);
    title(strcat(num2str(i + 2025),'年的暴雨数量预测值'));
    zlabel('暴雨次数 单位(次)');
    xlabel('经度 (°)');
    ylabel('纬度 (°)');
    
    colorbar;
    colormap(parula);
    view(0, 60);
    fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问3\预测结果\',strcat(num2str(i + 2025),'预测的结果'),'.jpg');
    saveas(gcf, fname);

end

%% 问4 根据人口和gdp、土地利用情况构建公式
%读取人口数组
popInfo = geotiffinfo('C:\Users\Administrator\Desktop\华为杯竞赛\附件数据集\数据集5\数据集5\数据集5\pop_output\pop1990\pop1990.tif');
popW = popInfo.CornerCoords.Lon(1:2);
popN = [2.07,56];
%terrain的图像经纬度范围
gdpInfo = geotiffinfo('C:\Users\Administrator\Desktop\华为杯竞赛\附件数据集\数据集6\数据集6\数据集6\GDP_output\gdp1990.tif');
gdpW = gdpInfo.CornerCoords.Lon(1:2);
gdpN = [0.88,57.8];

%% 读取人口数据,同样的数据大小不一致，池化处理
pop_data = zeros(1500, 800, 2015 - 1990 + 1);
pop_change = zeros(1,25);
for i = 1 : (2015 - 1990) + 1

    temper_str = dir(strcat('C:\Users\Administrator\Desktop\华为杯竞赛\附件数据集\数据集5\数据集5\数据集5\pop_output\','pop' ,num2str(1990 + i - 1)));
    temper_folder = string({temper_str.folder}');
    temper_file = string({temper_str.name}');
    temper_filepath = strcat(temper_folder, '\', temper_file);

    for j = 3 : length(temper_file)
        filepath_t = temper_filepath{j};
        if(contains(filepath_t, 'tif') && ~contains(filepath_t, 'aux') && ~contains(filepath_t, 'ovr'))
            pop_t = readTIF(filepath_t);
            disp(filepath_t);
            % temper_Matrix = temper_Matrix';
            temper_Matrix = pop_t;

            % pop_t(pop_t < 0) = NaN;
            pop_t(pop_t < 0) = 0;
            pop_change(i) = sum(pop_t(~isnan(pop_t)), 'all');

            pop_data(:,:,i) = MaxPool(temper_Matrix, 1500, 800);
        end
    end
end



%% 读取GDP文件,数据大小不一致，池化处理
GDP_data = zeros(1500, 800, 2015 - 1990 + 1);
GDP_change = zeros(1, 25);
GDP_change_rate = zeros(1, 24);
% imfinfo('C:\Users\Administrator\Desktop\华为杯竞赛\附件数据集\数据集6\数据集6\数据集6\GDP_output\gdp1990.tif');

temper_str = dir('C:\Users\Administrator\Desktop\华为杯竞赛\附件数据集\数据集6\数据集6\数据集6\GDP_output');
temper_folder = string({temper_str.folder}');
temper_file = string({temper_str.name}');
temper_filepath = strcat(temper_folder, '\', temper_file);
i = 1;
for j = 3 : length(temper_file)
    filepath_t = temper_filepath{j};
    
    if(contains(filepath_t, 'tif') && ~contains(filepath_t, 'xml') && ~contains(filepath_t, 'ovr') && ~contains(filepath_t, 'cpg') && ~contains(filepath_t, 'dbf'))
        disp(filepath_t);
        GDP_change_t = readTIF(filepath_t);
        temper_Matrix = GDP_change_t;
        
        % if(i ~= 26)
        %     GDP_change_t(GDP_change_t < 0) = NaN;
        % else
        %     GDP_change_t(GDP_change_t == 2.147483647000000e+09) = NaN;
        % end
         if(i ~= 26)
            GDP_change_t(GDP_change_t < 0) = 0;
        else
            GDP_change_t(GDP_change_t == 2.147483647000000e+09) = 0;
        end
        GDP_change(i) = sum(GDP_change_t(~isnan(GDP_change_t)), 'all');
        % temper_Matrix = temper_Matrix';
        GDP_data(:,:,i) = MaxPool(temper_Matrix, 1500, 800);
        i = i + 1;
    end
end
mesh(GDP_data(:,:,25)');

%% 土地利用率可以单指cropland
%提取一年的数据集进行训练
human_action = zeros(26 * 4184, 2);

 GDPW_step = (gdpW(2) - gdpW(1))/800;
 GDPN_step = (gdpN(2) - gdpN(1))/1500;


 popW_step = (popW(2) - popW(1))/800;
 popN_step = (popN(2) - popN(1))/1500;


for i = 1 : 26


    GDP_data_T = GDP_data(:,:,i)';
    pop_data_T = pop_data(:,:,i)';
    for j = 1 : 4184
        currentW = usedLandW(2) - ((127 - (row_usedLand(j))) * 0.5);
        currentN = usedLandN(2) - (col_usedLand(j) * 0.5);

        %提取GDP数据,和raincol一个比例
        GDPRow = round((currentW - gdpW(1)) / GDPW_step);
        GDPCol = round((gdpN(2) - currentN) / GDPN_step);

        human_action((i - 1) * 4184 + j, 1) = GDP_data_T(GDPRow, GDPCol);%第1个参数代表GDP

        %提取人口数据

        popRow = round((currentW - popW(1)) / popW_step);
        popCol = round((popN(2) - currentN) / popN_step);

        human_action((i - 1) * 4184 + j, 2) = pop_data_T(popRow, popCol);%第1个参数代表GDP

        if(i == 1 && j == 1)
            disp(currentW);
            disp(currentN);

            disp(GDPRow);
            disp(GDPCol);
            disp(popRow);
            disp(popCol);
        end

    end
end
%% GDP画图
figure(15);
GDP_data_T = GDP_data(:,:,25);
% GDP_data_T(GDP_data_T == 2.147483647000000e+09) = NaN;
GDP_data_T(GDP_data_T < 0) = NaN;

GDP_x_idx =  57.8 : -GDPN_step : 0.88 + GDPN_step;
GDP_y_idx =  52.010778892003614: GDPW_step : 1.544902172351913e+02 - GDPW_step;
mesh(GDP_y_idx, GDP_x_idx, GDP_data_T);
title('2015年全国GDP分布情况');
xlabel('纬度');
ylabel('经度');
zlabel('GDP值');
colorbar; % 显示颜色条以表示GDP值
view(0, 80);
% fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问4\','2015年全国GDP分布情况','.jpg');
% saveas(gcf, fname);

for i = 1 : 26
    if i > 1
        GDP_change_rate(i - 1) =  (GDP_change(i) - GDP_change(i - 1))/GDP_change(i - 1);
    end  
end

figure(17);
bar((1995:2015),GDP_change_rate(5:25),'FaceColor',[0.7 1 0.7])
title('1995年到2015年的GDP变化率');
xlabel('年份');
ylabel('幅度');
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问4\','1996年到2015年的GDP变化率','.jpg');
saveas(gcf, fname);

figure(18);
bar((1990:2015),GDP_change(1:26),'FaceColor', [1 1 0.8]);
title('1990年到2014年的GDP总量');
xlabel('年份');
ylabel('GDP总量');
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问4\','1990年到2015年的GDP总量','.jpg');
saveas(gcf, fname);
% 


%% 人口画图
figure(16);
pop_data_T = pop_data(:,:,25);

pop_data_T(pop_data_T < 0) = NaN;

pop_x_idx =  56 : -popN_step : 2.07 + popN_step;
pop_y_idx =  57.4892 : popW_step : 150.6537 - popW_step;
mesh(pop_y_idx,pop_x_idx,pop_data_T);
title('2015年全国人口分布情况');
xlabel('纬度');
ylabel('经度');
zlabel('人口密度');
colorbar; % 显示颜色条以表示GDP值
colormap(parula); % 可以尝试其他颜色映射如 'parula', 'hot', 'cool', 等等
view(0, 80);
% fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问4\','2015年全国人口分布情况','.jpg');
% saveas(gcf, fname);

pop_change_rate = zeros(1,25);
for i = 1 : 26
    if i > 1
        pop_change_rate(i - 1) =  (pop_change(i) - pop_change(i - 1))/pop_change(i - 1);
    end  
end

figure(18);
bar((1990:2015),pop_change,'FaceColor', [0.8 0.6 1]);
title('1990年到2015年的人口总数');
xlabel('年份');
ylabel('人口');
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问4\','1990年到2015年的人口总数','.jpg');
saveas(gcf, fname);

figure(19);
bar((1995:2015),pop_change_rate(5:25),'FaceColor', [0 0.298 0.133]);
title('1995年到2015年的人口变化数(相对上一年)');
xlabel('年份');
ylabel('人口变化率');
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问4\','1995年到2015年的人口变化率','.jpg');
saveas(gcf, fname);

%% 数据拼接
SVR_X = [SVMpredictionX(1:108784, 1:3), human_action, SVMpredictionX(1:108784, 4)];
%数据清洗
errIdx = zeros(1,26 * 4184);
idx_t = 1;
for i = 1 : (26 * 4184)
    if SVR_X(i, 1) == -111 || SVR_X(i, 2) == -32768 || SVR_X(i ,4) < 0 || SVR_X(i ,5) < 0 
        errIdx(idx_t) = i;
        idx_t = idx_t + 1;
    end
end
errIdx(errIdx==0) = [];
for j = length(errIdx) :-1: 1
    SVR_X([errIdx(j)],:) = [];
    SVR_X([errIdx(j)],:) = [];
end
%数据归一化
SVR_X_train = all2one(SVR_X(1: 3 * 4184, 1:5));
% SVR_X_train = SVR_X(1: 3 * 4184, 1:5);
SVR_Y_train = SVR_X(1: 3 * 4184, 6);
%测试集
SVR_X_Test = all2one(SVR_X(3 * 4184 : 4 * 4184, 1:5));
% SVR_X_Test = SVR_X(3 * 4184 : 4 * 4184, 1:5);
SVR_Y_Test = SVR_X(3 * 4184 : 4 * 4184, 6);
%% 相关性分析
% 计算相关系数矩阵
correlationMatrix = corrcoef(SVR_X_train);
% disp(correlationMatrix);
figure;
imagesc(correlationMatrix);
colorbar; % 添加颜色条
title('特征间的相关矩阵');
% xlabel('特征编号');

%% 训练
Model_q3 = fitrsvm(SVR_X_train, SVR_Y_train, 'KernelFunction', 'linear', 'BoxConstraint', 1);
% 定义超参数范围
% 定义超参数范围
% 定义超参数范围
C_values = logspace(-2, 2, 5); % 正则化参数 C
gamma_values = logspace(-2, 2, 5); % RBF 核的 gamma 参数

% 初始化最佳参数和最低误差
best_C = 0;
best_gamma = 0;
best_error = Inf;

% 网格搜索
for C = C_values
    for gamma = gamma_values
        fprintf('Trying C: %.4f, gamma: %.4f\n', C, gamma);
        mdl = fitrsvm(SVR_X_train, SVR_Y_train, 'KernelFunction', 'rbf', 'BoxConstraint', C, 'KernelScale', gamma);
        Y_pred = predict(mdl, SVR_X_Test);
        error = mean((Y_pred - SVR_Y_Test).^2); % 均方误差
        if error < best_error
            best_error = error;
            best_C = C;
            best_gamma = gamma;
            fprintf('New best: C: %.4f, gamma: %.4f, Error: %.4f\n', best_C, best_gamma, best_error);
        end
    end
end

% 打印最佳参数
fprintf('Best parameters found: C: %.4f, gamma: %.4f, Error: %.4f\n', best_C, best_gamma, best_error);

% 确保找到了有效的最佳参数
if best_C <= 0 || best_gamma <= 0
    error('最佳参数未找到或无效，请检查输入数据和参数范围。');
end

% 使用最佳参数训练最终模型
fprintf('Training final model with C: %.4f, gamma: %.4f\n', best_C, best_gamma);
final_mdl = fitrsvm(SVR_X_train, SVR_Y_train, 'KernelFunction', 'rbf', 'BoxConstraint', best_C, 'KernelScale', best_gamma);
predictedY = predict(final_mdl, SVR_X_Test);
% 计算模型的准确性
mse = mean((predictedY - SVR_Y_Test).^2); % 均方误差
rmse = sqrt(mse); % 均方根误差
rsquared = 1 - sum((SVR_Y_Test - predictedY).^2) / sum((SVR_Y_Test - mean(SVR_Y_Test)).^2); % 决定系数

% 显示结果
fprintf('均方误差 (MSE): %.4f\n', mse);
fprintf('均方根误差 (RMSE): %.4f\n', rmse);
fprintf('R平方: %.4f\n', rsquared);


plot(predictedY,'-o');
hold on;
plot(SVR_Y_Test,'-*');
xlabel('采样点');
ylabel('预测值');
set(gcf, 'WindowState', 'maximized');
title('在一年数据集下的预测表现');
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问4\','在一年数据集下的预测表现','.jpg');
saveas(gcf, fname);





figure;
residuals = SVR_Y_Test - predictedY;
scatter(predictedY, residuals, 'filled');
xlabel('预测值');
ylabel('残差');
title('残差图');
grid on;
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问4\','残差图','.jpg');
saveas(gcf, fname);
% 
figure;
scatter(SVR_Y_Test, predictedY, 'filled');
hold on;
plot(min(SVR_Y_Test):max(SVR_Y_Test), min(SVR_Y_Test):max(SVR_Y_Test), 'r--');
xlabel('实际值');
ylabel('预测值');
title('实际值 vs. 预测值');
legend('数据点', '理想线');
grid on;
fname = strcat('C:\Users\Administrator\Desktop\华为杯竞赛\问4\','实际值 vs. 预测值','.jpg');
saveas(gcf, fname);
hold off;