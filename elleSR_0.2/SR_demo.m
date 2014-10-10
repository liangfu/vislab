% SR_demo.m
%
% http://www.robots.ox.ac.uk/~elle/SRcode/
%
% (Just run this -- should be self-explanatory!)
%


% 1: Get a (synthetic) super-res problem:
[o,gtruth] = synthdata_demo; % Get struct "o" containing SR data, and "gtruth" ground truth image.
[biv,bih] = size(gtruth);
K = numel(o);


% 2: Look at the data:
suba = ceil(sqrt(3*K)/2);
subb = ceil(K/suba);
figure(1); clf;
for i = 1:K
    subplot(suba,subb,i); imagesc(o(i).im);
    axis image; colormap(gray); title(['low-res ' num2str(i)]); axis off;
end


% 3: Get matrices and set up options vector/params:
[W, Y, La, Lb, M] = makeW(biv,bih,o);
opts = zeros(1,18); % This is the "options" vector for the Netlab "scg" routine.
opts(1) = 1; % verbose
opts(2:3) = 1e-3; % convergeance criteria
opts(14) = 50; % number of iterations before automatic termination.
alp = 0.08;
nu = 0.04;


% 4: Find Average Image, Maximum Likelihood and Huber super-res image estimates:
[avim,msk,M] = getAvim(biv,bih,o);
im_ml = superres_ml(W,Y,La,Lb,avim,opts);
im_huber = superres_huber(W,Y,La,Lb,avim,alp,nu,opts);


% 5: Look at the various outputs, comparing them to the ground truth image:
figure(2); clf;
gap = 26;
subplot(2,2,1); imshow(avim((gap+1):end-gap,(gap+1):end-gap)+0.5);
axis image; colormap(gray); title('Average Image'); axis off;
subplot(2,2,2); imshow(im_ml((gap+1):end-gap,(gap+1):end-gap)+0.5);
axis image; colormap(gray); title('ML Image'); axis off;
subplot(2,2,3); imshow(im_huber((gap+1):end-gap,(gap+1):end-gap)+0.5);
axis image; colormap(gray); title('Huber Image'); axis off;
subplot(2,2,4); imshow(gtruth((gap+1):end-gap,(gap+1):end-gap)+0.5);
axis image; colormap(gray); title('Ground Truth Image'); axis off;


% 6: Be happy :-)
