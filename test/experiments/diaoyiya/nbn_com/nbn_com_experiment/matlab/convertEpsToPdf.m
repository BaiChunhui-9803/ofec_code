folderPath = 'H:\diaoyiya\ppt\ppt_20240924\PPT_20240926_nbn_tsp_analysis\experiments_show\nbn_eax_trait_figures/'; % 替换为实际的文件夹路径
epsFiles = dir(fullfile(folderPath, '*.eps'));

for i = 1:length(epsFiles)
    epsFileName = epsFiles(i).name;
    [~, epsBaseName, ~] = fileparts(epsFileName);
    pdfFileName = fullfile(folderPath, [epsBaseName '.pdf']);
    system(['epstopdf ', fullfile(folderPath, epsFileName), ' --outfile ', pdfFileName]);
end