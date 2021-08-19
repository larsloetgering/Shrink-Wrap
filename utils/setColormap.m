function r = setColormap

% cmap = [[0,0,0]; jet(20); [1,1,1]];
cmap = [[1,1,1]; jet(20)];

n = 100;
x = linspace(0,1,size(cmap, 1));
xq = linspace(0,1,n);

cmap2(:,1) = interp1(x,cmap(:,1),xq,'linear');
cmap2(:,2) = interp1(x,cmap(:,2),xq,'linear');
cmap2(:,3) = interp1(x,cmap(:,3),xq,'linear');

r = cmap2;
% cmap = [ [1,1,1];...
%          spectrumRGB(460);...
%          spectrumRGB(470);...
%          spectrumRGB(520);...
%          spectrumRGB(590);...
%          spectrumRGB(600);...
%          spectrumRGB(650);...
%          spectrumRGB(660);...
%          spectrumRGB(670)
%          ];
% n = 100;
% cmap2 = zeros(n, 3);
% 
% x = linspace(0,1,size(cmap, 1));
% xq = linspace(0,1,n);
% 
% cmap2(:,1) = interp1(x,cmap(:,1),xq,'linear');
% cmap2(:,2) = interp1(x,cmap(:,2),xq,'linear');
% cmap2(:,3) = interp1(x,cmap(:,3),xq,'linear');
% 
% % set return value
% r = cmap2;

end