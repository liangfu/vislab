function [o,gtruth] = synthdata_demo

% Create 5 low-res images using a built-in MATLAB image as ground truth.

gtruth = double(imread('westconcordorthophoto.png'))./255 - 0.5; % Load image and scale intensities to [-0.5,0.5]).

% westconcordorthophoto.png is a 366x364-pixel greyscale image that comes
% with Matlab's image processing toolbox. If you don't have it, try
% substituting any other greyscale image.

gam = 0.4;     % standard deviation of Gaussian point-spread function for the low-res images.
noise = 5/255; % standard deviation of noise on the low-res images.
zm = 2;        % zoom factor for the low-res images.
[biv,bih] = size(gtruth); % high-res image size.
K = 5;         % Number of low-res images to generate.

o = struct([]); % This will be my main super-res data structure.
for i = 1:K
    o(i).g = gam;    % standard deviation of Gaussian point-spread function.
    o(i).n = noise;  % standard deviation of noise.
    o(i).v = floor((biv*0.8)/zm);     % low-res image size (vertical).
    o(i).h = floor((bih*0.8)/zm);     % low-res image size (horizontal).
    o(i).la = gsamp(1,0.1^2,1);       % multiplicative photometric parameter (lambda_alpha)
    o(i).lb = gsamp(0,(10/255)^2,1);  % additive photometric parameter (lambda_beta)
    o(i).H = local_rand_homography(o(i).v,o(i).h,biv,bih,zm); % Invent a homography using function below.
    
    % Now generate a low-resolution image using the "gtruth" image and the
    % information we've just put into this structure:

    [o(i).im,o(i).noise,o(i).orig] = makeLR(gtruth,o(i));

end



return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = local_rand_homography(smv,smh,biv,bih,zm);
% This creates a random projective homography from a low-res image of size
% smv*smh to a high-res image of size biv*bih, with zoom factor zm.
%
T = diag([1/(smh-1),1/(smv-1),1]); % scaling
smshift = [eye(2),-0.5*[smh+1;smv+1];0,0,1]; % shift low-res image coords to be zero-centered.
bishift = [eye(2),0.5*[bih+1;biv+1];0,0,1];  % shift high-res image coords to be zero-centered.

P1 = [-0.5,0.5,0.5,-0.5;-0.5,-0.5,0.5,0.5;ones(1,4)];
P2 = [0.15*[rand(2,4)-0.5];zeros(1,4)]+P1;
stack = @(p,q) [zeros(1,3), -p(3)*q', p(2)*q'; p(3)*q', zeros(1,3), -p(1)*q'];
S = []; for i = 1:4, S = [S;stack(P1(:,i),P2(:,i))]; end;
H = bishift*diag([zm,zm,1])*inv(T)*reshape(null(S),3,3)'*T*smshift; 
H = H ./ H(9);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

