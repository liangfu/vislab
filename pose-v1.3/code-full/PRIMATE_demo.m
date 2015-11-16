function PRIMATE_demo
%clc; close all; clear;
globals;
name = 'PARSE';
% --------------------
% specify model parameters
% number of mixtures for 26 parts
K = 4*ones([1,14]);
% Tree structure for 26 parts: pa(i) is the parent of part i
% This structure is implicity assumed during data preparation
% (PRIMATE_data.m) and evaluation (PRIMATE_eval_pcp)
% pa = [0 1 2 3 4 5 6 3 8 9 10 11 12 13 2 15 16 17 18 15 20 21 22 23 24 25];
%pa = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
 pa = [0 1 2 3 4 5 3 3 3  4 10 11 13 14];
% Spatial resolution of HOG cell, interms of pixel width and hieght
% The PRIMATE dataset contains low-res people, so we use low-res parts
sbin = 4;
% --------------------
% Prepare training and testing images and part bounding boxes
% You will need to write custom *_data() functions for your own dataset
[pos neg test] = PRIMATE_data(name);pos
pos = point2box(pos,pa);pos
% --------------------
% training
try
  cls = [name '_final_' num2str(K')'];
  modeldata = load([cachedir cls]);model=modeldata.model;
catch
  model = trainmodel(name,pos,neg,K,pa,sbin);
  return
end

% --------------------
% testing phase 1
% human detection + pose estimation
suffix = num2str(K')';
model.thresh = min(model.thresh,-2);
boxes = testmodel(name,model,test,suffix);
if 0
det = PRIMATE_transback(boxes);
% --------------------
% evaluation 1: average precision of keypoints
% You will need to write your own APK evaluation code for your data structure
apk = eval_apk(det, test, 0.1);
% Average left with right and neck with top head
apk = (apk + apk([6 5 4 3 2 1 12 11 10 9 8 7 14 13]))/2;
% Change the order to: Head & Shoulder & Elbow & Wrist & Hip & Knee & Ankle
apk = apk([14 9 8 7 3 2 1]);
meanapk = mean(apk);
fprintf('mean APK = %.1f\n',meanapk*100);
fprintf('Keypoints & Head & Shou & Elbo & Wris & Hip  & Knee & Ankle\n');
fprintf('APK       '); fprintf('& %.1f ',apk*100); fprintf('\n');
end
% --------------------
% testing phase 2
% pose estimation given ground truth human box
model.thresh = min(model.thresh,-2);
boxes_gtbox = testmodel_gtbox(name,model,test,suffix);
if 0
det_gtbox = PRIMATE_transback(boxes_gtbox);
% --------------------
% evaluation 2: percentage of correct keypoints
% You will need to write your own PCK evaluation code for your data structure
pck = eval_pck(det_gtbox, test, 0.1);
% Average left with right and neck with top head
pck = (pck + pck([6 5 4 3 2 1 12 11 10 9 8 7 14 13]))/2;
% Change the order to: Head & Shoulder & Elbow & Wrist & Hip & Knee & Ankle
pck = pck([14 9 8 7 3 2 1]);
meanpck = mean(pck);
fprintf('mean PCK = %.1f\n',meanpck*100); 
fprintf('Keypoints & Head & Shou & Elbo & Wris & Hip  & Knee & Ankle\n');
fprintf('PCK       '); fprintf('& %.1f ',pck*100); fprintf('\n');
end
% --------------------
% visualization
% figure(1);
% visualizemodel(model);
% figure(2);
% visualizeskeleton(model);
demoimid = 7;
im = imread(test(demoimid).im);
%colorset = {'g','g','y','r','r','r','r','y','y','y','m', ...
%'m','m','m','y','b','b','b','b','y','y','y','c','c','c','c'};
colorset={};for ii=1:6,colorset{ii}='r';end;
for ii=7:12,colorset{ii}='g';end;for ii=13:14,colorset{ii}='k';end
box = boxes{demoimid};
% show all detections
figure(3);
subplot(1,2,1); showboxes(im,box,colorset);
subplot(1,2,2); showskeletons(im,box,colorset,model.pa);
% show best detection overlap with ground truth box
box = boxes_gtbox{demoimid};
figure(4);
subplot(1,2,1); showboxes(im,box,colorset);
subplot(1,2,2); showskeletons(im,box,colorset,model.pa);