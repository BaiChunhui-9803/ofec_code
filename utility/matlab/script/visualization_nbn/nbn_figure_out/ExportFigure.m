function []= ExportFigure(fid,filename, filetype)
   filenameStr=sprintf('%s',filename);
   filetypeStr=sprintf('%s',filetype);
   fidStr=sprintf('%d',fid);
   fh = figure(fid);
   hold on;
   fh.WindowState = 'maximized';
%   fh.Position = [0,0,1080,1080];
   axis square;
   view([180 0]);
   filepath=[filenameStr,'_',fidStr,'_front.',filetypeStr];
   exportgraphics(gcf,filepath);
   view([-37.5000 30]);
   filepath=[filenameStr,'_',fidStr,'_origin.',filetypeStr];
   exportgraphics(gcf,filepath);
   view(2);
   filepath=[filenameStr,'_',fidStr,'_top.',filetypeStr];
   exportgraphics(gcf,filepath);
   hold off;
end 