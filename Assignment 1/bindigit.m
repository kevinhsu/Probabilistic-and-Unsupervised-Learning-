% Assignment 3, Unsupervised Learning, UCL 2003
% Zoubin Ghahramani
% Matlab comments start with %

% load the data set

load binarydigits.txt -ascii;
Y=binarydigits;
[N D]=size(Y);

% this is how you display one image, e.g. the 4th image:
%
% y=Y(4,:);
% colormap gray;
% imagesc(reshape(y',8,8)')

% you can also reshape by hand as follows (slower, but useful for
% non-Matlab implementations)
% for i=1:8,
%  for j=1:8,
%    x(i,j)=y((i-1)*8+j);
%  end;
% end;

% now we will display the whole data set:
colormap gray;
for n=1:N, 
  subplot(10,10,n);
  imagesc(reshape(Y(n,:)',8,8)'); 
  axis off;
end;

