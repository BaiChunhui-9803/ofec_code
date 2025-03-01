function []= setExportFigureType(filename,veiwName,offset)
   filenameStr=sprintf('%s',filename);
 %  filetypeStr=sprintf('%s',filetype);
   veiwNameStr=sprintf('%s',veiwName);
   ax = gca;
   ax.OuterPosition = [0 0 1 1];
   outerpos = ax.OuterPosition; % 获取外部框位置
   ti = ax.TightInset; % 获取内容框位置
   left = outerpos(1) + ti(1); % 把Position的left值设为左边margin的值
   bottom = outerpos(2) + ti(2) ; % 把Position的bottom值设为右边margin的值
   ax_width = outerpos(3) - ti(1) - ti(3); % 设置对应的宽度
   ax_height = outerpos(4) - ti(2) - ti(4); % 设置对应的高度
   ax.Position = [left bottom ax_width-offset ax_height-offset]; % 可以微调一下，以保证边缘没有被剪裁掉。
   fig = gcf;
   fig.PaperPositionMode = 'auto';
   fig_pos = fig.PaperPosition;
   fig.PaperSize = [fig_pos(3) fig_pos(4)];
%   filepath=[filenameStr,'-',veiwNameStr,'.',filetypeStr];
%   saveas(gcf,filepath,fileSuffix);
%   filepath=[filenameStr,'-',veiwNameStr,'.pdf'];
%   exportgraphics(gcf,filepath,'ContentType','vector');
 %     filepath=[filenameStr,'-',veiwNameStr,'.eps'];
 %  exportgraphics(gcf,filepath,'ContentType','vector');
   filepath=[filenameStr,'-',veiwNameStr,'.png'];
   exportgraphics(gcf,filepath);
   %disp(outerpos);
   %disp(ti);
   ax.OuterPosition = [0 0 1 1];
end