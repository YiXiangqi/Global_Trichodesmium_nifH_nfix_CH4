function rgbMatrix = redToWhiteColorbarCustom(n)
    % redToWhiteColorbarCustom 生成从红色到白色的自定义渐变色彩条
    % 输入 n：需要的颜色数量
    % 返回 n 行 3 列的 RGB 矩阵
    
    if n < 2 || n > 256
        error('n must be between 2 and 256');
    end
    
    % 样本颜色的 RGB 数据（按行排列）
    sampleColors = [
        92     0    15;
        141    15    20;
        180    20    25;
        212    33    32;
        239    60    44;
        250    94    64;
        252   127    95;
        252   161   132;
        253   195   172;
        254   224   211
    ];
    
    % 归一化样本颜色到 [0, 1]
    sampleColors = sampleColors / 255;
    
    % 样本位置（线性分布）
    samplePositions = linspace(0, 1, size(sampleColors, 1));
    
    % 目标位置（根据 n 平均分布）
    targetPositions = linspace(0, 1, n);
    
    % 插值每个通道的颜色值
    r = interp1(samplePositions, sampleColors(:, 1), targetPositions);
    g = interp1(samplePositions, sampleColors(:, 2), targetPositions);
    b = interp1(samplePositions, sampleColors(:, 3), targetPositions);
    
    % 合成 RGB 矩阵
    rgbMatrix = [r', g', b'];
end
