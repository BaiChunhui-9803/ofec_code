function []= ExportThreeFigures(filename)
   fh= gcf;
   fh.WindowState = 'maximized';
   view([-37.5000 30]);
   setExportFigureType(filename,'origin',0.1);
   view(2);
   setExportFigureType(filename,'top',0.1);
      view([180 0]);
   setExportFigureType(filename,'front',0.1);
%      view([-37.5000 30]);
%   setExportFigureType(filename,'origin',filetype,fileSuffix,0.1);
end 