function []= ExportFigures(figureIds,filename, filetype)
   filenameStr=sprintf('%s',filename);
   filetypeStr=sprintf('%s',filetype);
for fid = figureIds
   fidStr=sprintf('%d',fid);
   fh = figure(fid);
%fh.WindowState = 'maximized';
   view([180 0]);
   filepath=[filenameStr,'_',fidStr,'_front.',filetypeStr];
exportgraphics(gcf,filepath);
   view([-37.5000 30]);
   filepath=[filenameStr,'_',fidStr,'_origin.',filetypeStr];
   exportgraphics(gcf,filepath);
   view(2);
   filepath=[filenameStr,'_',fidStr,'_top.',filetypeStr];
   exportgraphics(gcf,filepath);
end
end 