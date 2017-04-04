function [so, spatialRFonedim, magnitude1STD] = buildSpatialRF(rfDiameter)

% rfDiameter = rfDiameter / (patchSize(2) / inCol);

extent = 5;    % ratio between sampling size and spatial RF standard dev
r = 0.75;        % radius ratio between center and surround for DoG
k = 1.032 * r;   % scaling of magnitude of surround


% number bipolar cells out to the extent of the spatial RF
pts = -extent*rfDiameter+1 : extent*rfDiameter;

% Specify centers in um, offset even rows for hexagonal packing
% ic = centerX(ii) - (mod(jj, 2) - 0.5) * 2* rfDiameter + 0*3*centerNoise*(2*rand(1,1)-1);
% jc = centerY(jj) + 0*3*centerNoise*(2*rand(1,1)-1);
ic = rfDiameter/2; jc = rfDiameter/2;

% Add some noise to deviate from circularity
% (unitless: Q = (1/d^2)*[1 0; 0 1] yields circular SD with r = d
d1 = 1; d2 =  0*10*0.0675*(rand(1,1)-0.5);      % 0.0675*randn(1,1);
Q = (1/rfDiameter^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);

% Calculate values for input to DoG function in an efficient way
[i2, j2] = meshgrid(ic+pts, jc+pts); % um
i = i2(:); j = j2(:);                % um

IJ = bsxfun(@minus,[i j],[ic jc]); % um
QIJ = Q*IJ'; rQIJ = r*Q*IJ';       % unitless
%  icrm = repmat([ic jc],length(i),1);

% (-0.5*(x-c)*Q*(x-c)'): unitless
p1 = prod([IJ(:,1) QIJ(1,:)'],2)+prod([IJ(:,2) QIJ(2,:)'],2);
p2 = prod([IJ(:,1) rQIJ(1,:)'],2)+prod([IJ(:,2) rQIJ(2,:)'],2);

% DoG calculation
% conditional intensity, related by Poisson firing to spikes/sec
so_center = reshape(exp(-0.5*p1), size(i2));
so_surround = reshape(k*exp(-0.5*p2), size(i2));
so = so_center - so_surround;


%% Vectors (1D) instead of matrices (2D)
% conditional intensity, related by Poisson firing to spikes/sec
sx_cent = exp(-0.5*Q(1,1)*(0+pts).^2);
sy_cent = exp(-0.5*Q(2,2)*(0+pts).^2);
sx_surr = sqrt(k)*exp(-0.5*Q(1,1)*r*(0+pts).^2);
sy_surr = sqrt(k)*exp(-0.5*Q(2,2)*r*(0+pts).^2);

% Store calculated parameters, units of conditional intensity
cellCenterLocations = [ic jc];% - [centerCorrectX centerCorrectY]; % um
spatialRFArray = so;
sRFcenter = so_center;
sRFsurround = so_surround;

% units of conditional intensity
spatialRFonedim = [(sx_cent - sx_surr); (sy_cent - sy_surr)];

%% Measure contour at 1 SD
% Choose random orientation of asymmetrical RF to measure 1 SD
% magnitude
xv = [1 0];%rand(1,2);
xvn = rfDiameter*xv;% * xv./norm(xv);
x1 = xvn(1); y1 = xvn(2);

% Do some calculations to make plots where RFs are filled in
% Measure magnitude at 1 SD from center
magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
% Find components of RF over that magnitude and store
spatialRFFill  = find(abs(so_center)>magnitude1STD);
rfDiaMagnitude = magnitude1STD;
