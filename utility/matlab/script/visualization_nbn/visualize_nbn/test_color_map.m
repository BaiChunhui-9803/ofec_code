x=rand(20,1); %data to be plotted
ran=range(x); %finding range of data
min_val=min(x);%finding maximum value of data
max_val=max(x);%finding minimum value of data
y=floor(((x-min_val)/ran)*199)+1; 
col=zeros(20,3);
%p  = colormap(jet(200));
map = jet(200);
for i=1:20
  a=y(i);
  col(i,:)= map(a,:);
  stem3(i,i,x(i),'Color',col(i,:))
  hold on
end