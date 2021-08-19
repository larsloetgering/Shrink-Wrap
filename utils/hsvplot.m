function varargout = hsvplot(u, varargin)
% hsvplot(u) generates hue-brightness plot of two dimensional input u
% last change: 3rd March 2018

p = inputParser;
p.addParameter('intensityScale', [ ])
p.addParameter('colorbar', false);
p.addParameter('scalebar', 0);
p.parse( varargin{:} );

% normalize birghtness (value) to range [0, 1]
u = gather(u);
r = abs(u);
r = r / ( max(r(:)) + eps );
if ~isempty(p.Results.intensityScale)
    r = posit(r - p.Results.intensityScale(1));
    r = r/(p.Results.intensityScale(2) - p.Results.intensityScale(1));
end

% normalize angle
phi = angle( u );
phi = ( phi + pi )/( 2 * pi );

% normalization of phase saturation

B = zeros(size(u,1), size(u,2), 3);         %Declare RGB array
B(:,:,1) = phi;
B(:,:,2) = 1;
B(:,:,3) = r;
A = hsv2rgb(B);

switch p.Results.colorbar
    
    case 'both'
        % draws both hue and brightness bars (left and bottom)
        
        k = ceil(size(u,1)/15);
        
        % phase bar
        m = size(u,1);
        D = zeros(m, 1, 3);
        D(:,:,1) = (0:1/m:1-1/m)';
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        C(:,1:k,:) = repmat( C(:,1,:),1, k );
        A=[C, A];
        
        % intensity bar
        n = size(A, 2);
        E = zeros(1, n, 3);
        E(:,:,1) = 1;
        E(:,:,2) = 0;
        E(:,:,3) = (0:1/n:1-1/n);
        F = hsv2rgb(E);
        F = repmat( F, k, 1 );
        A=[A; F];
        
    case 'North'
        
        m = size(u,1);
        k = ceil(m/15);
        
        % phase bar
        D = zeros(1, m, 3);
        D(:,:,1) = (0:1/m:1-1/m);
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        C = repmat( C, k, 1 );
        A(1:k, :, :) = C;
        
    case 'East'
        
        m = size(u,1);
        k = ceil(m/15);
        
        % phase bar
        D = zeros(m, 1, 3);
        D(:,:,1) = (0:1/m:1-1/m)';
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        C = repmat( C, 1, k );
        A(:, m-k+1:m, :) = C;
        
    case 'South'
        
        m = size(u,1);
        k = ceil(m/15);
        
        % phase bar
        D = zeros(1, m, 3);
        D(:,:,1) = (0:1/m:1-1/m);
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        C = repmat( C, k, 1 );
        A(m-k+1:m, :, :) = C;
        
    case 'West'
        
        m = size(u,1);
        k = ceil(m/15);
        
        % phase bar
        D = zeros(m, 1, 3);
        D(:,:,1) = (0:1/m:1-1/m)';
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        C = repmat( C, 1, k );
        A(:, 1:k, :) = C;
        
        
        
end

if p.Results.scalebar > 0
    
    m = round(size(A, 1)/20);
    n = round(size(A, 2)/20);
    scaleBar = ones( m + 1, p.Results.scalebar + 1, 3 ,'like', A);
    scaleBar(:,:,1) = 0.1;
    A(m:2*m,(3*n):(3*n+p.Results.scalebar),:) = hsv2rgb(scaleBar);
    
end

imagesc(A); axis image off
set(gcf, 'Color', 'w');

switch nargout
    case 1
        varargout{1}=uint8(A*255);
    otherwise
end

return