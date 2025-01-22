
filedir = '//172.29.71.217/e/DiaoYiya/a_final_NBN_result/neutrality/';
filename = 'ic_compare10';
filepath=[filedir,filename,'.txt'];
fileID = fopen(filepath,'r');

tline = fgetl(fileID);
formatSpec = '%e';
mat = fscanf(fileID,formatSpec,[9,99]);
x= mat(2,:);
hmax = mat(3,:);
eps_s = mat(4,:);
eps_max = mat(5,:);
eps_ratio = mat(6,:);
mo = mat(7,:);
nbn = mat(8,:);



trace1 = struct(...
  'x', x, ...
  'y', hmax, ...
  'name', 'h.max', ...
  'type', 'scatter');

trace2 = struct(...
  'x', x, ...
  'y', eps_s, ...
  'name', 'eps.s', ...
  'yaxis', 'y2', ...
  'type', 'scatter');

trace3 = struct(...
  'x', x, ...
  'y', eps_max, ...
  'name', 'eps.max', ...
  'yaxis', 'y3', ...
  'type', 'scatter');

trace5 = struct(...
  'x', x, ...
  'y', mo, ...
  'name', 'mo', ...
  'yaxis', 'y5', ...
  'type', 'scatter');



trace6 = struct(...
  'x', x, ...
  'y', eps_ratio, ...
  'name', 'eps.ratio', ...
  'yaxis', 'y6', ...
  'type', 'scatter');



trace4 = struct(...
  'x', x, ...
  'y', nbn, ...
  'name', 'nbn', ...
  'yaxis', 'y4', ...
  'type', 'scatter');

data = {trace1, trace2, trace3,trace5,trace6, trace4};

layout = struct(...
    'title', 'neutrality-dim10', ...
    'width', 800, ...
    'xaxis', struct('domain', [0.3, 0.7]), ...
    'yaxis', struct(...
      'titlefont', struct('color', '#0072BD'), ...
      'tickfont', struct('color', '#0072BD')), ...
    'yaxis2', struct(...
      'titlefont', struct('color', '#D95319'), ...
      'tickfont', struct('color', '#D95319'), ...
      'anchor', 'free', ...
      'overlaying', 'y', ...
      'side', 'left', ...
      'position', 0.15), ...
    'yaxis3', struct(...
      'titlefont', struct('color', '#EDB120'), ...
      'tickfont', struct('color', '#EDB120'), ...
      'anchor', 'x', ...
      'overlaying', 'y', ...
      'side', 'right'), ...
          'yaxis5', struct(...
      'titlefont', struct('color', '#7E2F8E'), ...
      'tickfont', struct('color', '#7E2F8E'), ...
      'anchor', 'x', ...
      'overlaying', 'y', ...
      'side', 'right'), ...
          'yaxis6', struct(...
      'titlefont', struct('color', '#77AC30'), ...
      'tickfont', struct('color', '#77AC30'), ...
      'anchor', 'x', ...
      'overlaying', 'y', ...
      'side', 'right'), ...
    'yaxis4', struct(...
      'titlefont', struct('color', '#4DBEEE'), ...
      'tickfont', struct('color', '#4DBEEE'), ...
      'anchor', 'free', ...
      'overlaying', 'y', ...
      'side', 'right', ...
      'position', 0.85));

plotly(data, struct('layout', layout));
