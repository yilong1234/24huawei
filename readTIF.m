function [outputMatrix] = readTIF(filepath)
    %filepath='C:\Users\Administrator\Desktop\华为杯竞赛\中国大陆0.5°土地利用和覆盖变化数据集(1900-2019年)\数据实体\cropland-1900.tif';                                         %%图像名称与路径
    Info=imfinfo(filepath);                                      %%获取图片信息并判断是否为tif
    
    tif='tif';
    format=Info.Format;
    if  (strcmp(format ,tif)==0)
        disp('载入的不是tif图像，请确认载入的数据');                %%确保载入的图像是tiff图像
    end
    
    Slice=size(Info,1);                                          %%获取图片z向帧数
    Width=Info.Width;
    Height=Info.Height;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Image=zeros(Height,Width,Slice);
    
    for i=1:Slice
        Image(:,:,i)=imread(filepath,i);                         %%一层一层的读入图像
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for i=1:Slice
    %     J=uint8(Image(:,:,i));                                   %%一层一层写出图像
    %     %%imwrite(J,[num2str(i,'%4d'),'.tif'],'WriteMode','Append');
    %     imwrite(J,[num2str(i,'%04d'),'.tif']);
    % end
    
    % mesh(Image);
    % if(year >= 1990 && year <= 2020)
    %     pcolor(Image);
    %     % 设置坐标轴刻度
    %     axis ij;
    %     axis tight;
    %     % 移除网格线
    %     clim([0 1]); % 设置颜色条的最小值和最大值
    %     set(gca, 'GridLineStyle', 'none'); % 关闭网格线
    % 
    %     % 设置标题和其他属性
    %     titleValue = strcat('土地利用情况：',type);
    %     title(titleValue);
    %     xlabelStr = strcat('当前年份为', num2str(year));
    %     xlabel(xlabelStr);
    %
    %     % 显示颜色条
    %     colormap(jet); % 可以选择不同的颜色图
    %     colorbar;
    %
    %       outputMatrix = Image;
    % end
    outputMatrix = Image;
end
