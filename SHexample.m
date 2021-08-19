%% shrink-wrap 2D example

% description: this is a toy code implementation of shrink wrap
% original work by 
% Marchesini, Stefano, et al. "X-ray image reconstruction from a diffraction 
% pattern alone." Physical Review B 68.14 (2003): 140101.
% 
% code written by: Lars Loetgering
% last update: 19th August 2021
% 
% notes: 
% - here a near-field propagator (aspw) is used. Change propagator to fft
% to simulate Fraunhofer diffraction geometry

% configurations
clear
close all

set(0,'DefaultFigureWindowStyle','normal')
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
addpath(genpath('utils'))
set(0,'defaultfigurecolor',[1,1,1])

%% step 1) create world-shaped, halftone object ( = object)
% read object
load('world.mat');

% create checkerboard for dots
N = 2^9;
im1 = zeros(N,'single');
k = 30;
im1(1:30:end,:) = 1;
im1(k/2:k:end,:) = 1;

im2 = zeros(N,'single');
im2(:,1:k:end) = 1;
im2(:,k/2:k:end,:) = 1;

im3 = im1.*im2;
idx = find(abs(im3(:))==1);

% physical parameters
lambda = 632.8e-9;
dx = 10e-6;
L = N * dx;
x = (-N/2:N/2-1)*dx;
[X,Y] = meshgrid(x);
z = 10e-2;   % propagation distance for ASP propagator

im = imresize(single(world), N/2*[1,1]);

im = diff_mask(im, N*[1,1]);
im = real(center(im));
object = single(im3.*im);

% shuffle some of the phases of the dot-shaped regions
rng(1)
for k = 1:length(idx)
    if rand(1,1) > 0.7
        object(idx(k)) = object(idx(k)) * exp(1i * 2*pi*rand(1,1));
    end
end

object = convolve2(object, fspecial('disk',5),'same');
object = object / max(abs(object(:)));
object = (abs(object) >= 0.5) .* exp(1i*angle(object));
disp(['nnz in object: ', num2str(sum(abs(object(:))))])

% show image
figure(1)
n = 256;
subplot(1,2,1)
hsvplot(object((N/2-n+1:N/2+n), (N/2-n+1:N/2+n)),'colorbar','North')
set(gcf, 'Color', 'w');
title('object')

subplot(1,2,2)
hsvplot(conj(object((N/2-n+1:N/2+n), (N/2-n+1:N/2+n))),'colorbar','North')
set(gcf, 'Color', 'w');
title('conj(object)')
% note: both the object and it's complex conjugate are possible solutions
% when a far field propagator is used (to this end, replace aspw with fft)

%% step 2) simulate diffraction data

bitDepth = 2^14;
diffPattern = abs( aspw(object, z, lambda, L) ).^2;
diffPattern = diffPattern / max(diffPattern(:)) * bitDepth;

% add poisson noise
diffPattern(:) = poissrnd(diffPattern(:));
cmap = setColormap;
figure(2)
imagesc(x,x,log10(diffPattern+1))
axis image off; colormap(cmap); colorbar
set(gcf, 'Color', 'w');
title(['diffraction pattern @ z=',num2str(z)])

% convert to modulus
diffPattern = sqrt(diffPattern);

%% step 3) initialization and preallocation
% initial support estimate

S = circ(X, Y, 300*dx);

figure(3)
imshowpair(S, abs(object))
axis image; title('initial support and object')

% initial object estimate

g = aspw(diffPattern, -z, lambda, L);
figure(4)
hsvplot(g)

%% preallocate transfer function for angular spectrum propagator

gpuBool = true;
[~,H] = aspw(g, z, lambda, L);

try
    if gpuBool
        g = gpuArray(g);
        H = gpuArray(H);
    end
catch
    error('set gpuBool = false if your computer has no gpu available.')
end

%% error

errorMetric = [ ];

%% step 4) shrink wrap algorithm

noi = 100000;                   % number of iterations
swFrequency = 5;                % determines frequency of support re-estimation
figureUpdateFrequency = 1000;   % determines frequency of plots
for k = 1:noi
    
    % 1) propagate estimate, g, to detector plane
    G = ifft2c(fft2c(g) .* H);
    
    % calculate error between estimate and measurement
    errorMetric = [errorMetric, sum( (abs(G(:)) - diffPattern(:)).^2 )];
    
    % 2) replace modulus by measured diffraction , keep phase
    G = diffPattern .* exp(1i * angle(G));
    
    % 3) backpropagate to object plane
    gUpdate = ifft2c(fft2c(G) .* conj(H));
    
    % 4) obtain improved mask (shrink-wrap operation)
    if mod(k, swFrequency) == 0 || k == 1
        
        S = abs(gUpdate); % make real-valued
        %         h = fspecial('gaussian', 20, 2);
        h = fspecial('disk', 5);
        S = convolve2(S, h, 'same');
        S = S/max(S(:)); % make in range [0,1]
        
        % shrink
        thresh = 0.15;
        S = S > thresh;
        
    end
    
    % 4) relaxed averaged succesive reflection (RAAR) algorithm including
    % support constraint
    gp = 2*gUpdate - g;
    gp = 2*S.*gp - gp;
    beth = 0.5;
    g = beth * (gp+g)/2 + (1-beth)*gUpdate;
    
    gimmel = 1e-2;
    g = gimmel*abs(g) + (1-gimmel)*g;
    
    % plot estimtae
    if mod(k, figureUpdateFrequency) == 0 || k == 1
        
        temp = phaseSynchronization(object, center(g));
        
        figure(5); imagesc(x(N/2-n+1:N/2+n), x(N/2-n+1:N/2+n), S((N/2-n+1:N/2+n), (N/2-n+1:N/2+n)))
        axis image off, colormap gray
        title({'support estimate',['nonzero elements:',num2str(sum(S(:)))]})
        
        figure(6)
        hsvxplot(g((N/2-n+1:N/2+n), (N/2-n+1:N/2+n)),'intensityScale',[0 0.7])
        title(['iteration: ', num2str(length(errorMetric))])
        
        figure(7)
        loglog(errorMetric(2:end),'o')
        xlabel('iteration'), ylabel('error')
        grid on
        drawnow
        
    end
    
end