function ploterp(X,X2,offset,times,color)
% ploterp(X,X2,offset,time,color) X and X2 should be trials by
% samples. offset draws a line at that time. times is a vector that
% tells the times of each sample for plotting. color is a 3-vector
% indicating rgb values for the line and SEM area.

if nargin<4 | isempty(times), times=1:size(X,2); end
if nargin<5 | isempty(color), 
  linecolor = [0 0 1]; backcolor=[0.8 0.8 0.95]; 
else
  linecolor = color; backcolor=[0.8 0.8 0.8]+color/5; 
end

[N,L] = size(X);

if nargin<2 | isempty(X2)
  m = mean(X);
  se = std(X)/sqrt(size(X,1));
else
  m = mean(X)-mean(X2);
  se = sqrt((std(X)/sqrt(size(X,1))).^2+(std(X2)/sqrt(size(X,1))).^2);
end

plot(times([1 L]), [0 0],'k'); hold on;

fill([times(1:L)  times(L:-1:1)],[m+2*se m(end:-1:1)-2*se(end:-1:1)], ...
    backcolor,'LineStyle','n')
axis tight; 
ax=axis;
plot(times,m,'color',linecolor);    
hold off

if nargin>2 & ~isempty(offset) 
  hold on;
  plot([offset offset],ax(3:4),'k')
  axis(ax);  
  hold off
end