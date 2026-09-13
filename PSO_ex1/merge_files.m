function merge_files(outputFileName)
% MERGE_M_FILES 合并当前工作文件夹下所有 .m 文件为一个文本文件
%   使用方法：
%   1. merge_m_files          -> 默认生成 'combined_code.txt'
%   2. merge_m_files('my.txt') -> 自定义生成文件名

    % 如果用户没有输入文件名，使用默认名称
    if nargin < 1
        outputFileName = 'combined_code.txt';
    end

    % 关键点：获取当前你“正在浏览/工作”的文件夹路径（pwd）
    % 这样无论函数存在哪里，它操作的都是你当前眼前的文件夹
    currentFolder = pwd;
    mFiles = dir(fullfile(currentFolder, '*.m'));

    if isempty(mFiles)
        disp('当前文件夹下没有找到任何 .m 文件。');
        return;
    end

    % 在当前工作目录下创建输出文件
    fullOutputPath = fullfile(currentFolder, outputFileName);
    fidOut = fopen(fullOutputPath, 'w');

    if fidOut == -1
        error('无法创建或打开输出文件 %s', fullOutputPath);
    end

    % 循环遍历当前文件夹下的所有 .m 文件
    for k = 1:length(mFiles)
        currentFile = mFiles(k).name;
        
        % 如果输出文件名正好也是 .m 或者是重名的 txt，跳过它防止死循环
        if strcmp(currentFile, outputFileName)
            continue;
        end
        
        fprintf('正在合并: %s...\n', currentFile);
        
        % 写入文件头信息
        fprintf(fidOut, '==================================================\n');
        fprintf(fidOut, '%% 文件名: %s\n', currentFile);
        fprintf(fidOut, '==================================================\n\n');
        
        % 读取并写入内容
        try
            fileContent = fileread(fullfile(currentFolder, currentFile));
            fprintf(fidOut, '%s\n', fileContent);
        catch
            fprintf('读取失败: %s\n', currentFile);
        end
        
        % 每个文件结束后换两行
        fprintf(fidOut, '\n\n');
    end

    % 关闭输出文件
    fclose(fidOut);
    fprintf('完成！已在当前目录下生成合并文件: %s\n', outputFileName);
end