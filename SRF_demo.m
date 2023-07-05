clear;
close all
addpath(genpath(pwd))
%% load dataset
% dataset#1 to dataset#6
dataset = 'dataset#1';
Load_dataset % For other datasets, we recommend a similar pre-processing as in "Load_dataset"
fprintf(['\n Data loading is completed...... ' '\n'])
%% Parameter setting
% With different parameter settings, the results will be a little different
% Ns: the number of superpxiels,  Ns = 5000 is recommended.
% eta: 0.1 <= eta <= 0.9 is recommended.
%--------------------Available datasets and Parameters---------------------
%    Dataset       |                eta               |
%   dataset#1      |                0.5               |
%   dataset#2      |                0.6               |
%   dataset#3      |                0.5               |
%   dataset#4      |                0.5               |
%   dataset#5      |                0.3               |
%   dataset#6      |                0.3               |
%--------------------------------------------------------------------------
opt.Ns = 5000;
opt.lambda = 0.1;
opt.beta =1;
opt.eta = 0.5;
opt.alfa = 0.05; % MRF segmentation parameter
opt.fuse_type = 'dot_fixed_Gradient';% 
if strcmp(dataset,'dataset#1') == 1 || strcmp(dataset,'dataset#3') == 1 || strcmp(dataset,'dataset#4') == 1
    opt.eta = 0.5;
elseif strcmp(dataset,'dataset#2') == 1
    opt.eta = 0.6;
    opt.alfa = 0.2;
elseif strcmp(dataset,'dataset#5') == 1
    opt.eta = 0.3;
elseif strcmp(dataset,'dataset#6') == 1
    opt.eta = 0.3;
    opt.fuse_type = 'dot_fixed';
end
%%
fprintf(['\n SRF is running...... ' '\n'])
time = clock;
[regression_t1,regression_t2,DI_t1,DI_t2,CM] = SRF_main(image_t1,image_t2,opt);
fprintf(['\n' '====================================================================== ' '\n'])
fprintf('\n');fprintf('The total computational time of SRF is %i \n' ,etime(clock,time));
fprintf(['\n' '====================================================================== ' '\n'])
%% Displaying results
fprintf(['\n Displaying the results...... ' '\n'])
[tp,fp,tn,fn,fplv,fnlv,~,~,pcc,kappa,imw]=performance(CM,Ref_gt);
F1 = 2*tp/(2*tp+fp+fn);
figure;
subplot(331);imshow(image_t1);title('original imgt1')
subplot(332);imshow(image_t2);title('original imgt2')
subplot(333);imshow(Ref_gt);title('ground truth')
subplot(334);imshow(uint8(regression_t1));title('regression imgt1')
subplot(335);imshow(uint8(regression_t2));title('regression imgt2')
subplot(336);imshow(CM);title('change map')
subplot(337);imshow(remove_outlier(DI_t1),[]);title('backward DI')
subplot(338);imshow(remove_outlier(DI_t2),[]);title('forward DI')
subplot(339);imshow(CMplotRGB(CM,Ref_gt));title('change map')

if strcmp(dataset,'dataset#4') == 1 | strcmp(dataset,'dataset#5') == 1 % for SAR image in dataset#4 and #5
    subplot(334);imshow(uint8(exp(regression_t1*(log(256)-log(1)))-1));title('regression imgt1')
elseif strcmp(dataset,'dataset#6') == 1
    subplot(335);imshow(uint8(exp(regression_t2*(log(255)-log(6))+log(6))-1));title('regression imgt2')
end
[tp,fp,tn,fn,fplv,fnlv,~,~,OA,kappa]=performance(CM,Ref_gt);
F1 = 2*tp/(2*tp+fp+fn);
result = 'SRF: OA is %4.3f; Kc is %4.3f; F1 is %4.3f \n';
fprintf(result,OA,kappa,F1)
%%
if F1 < 0.3
    fprintf('\n');disp('Please select the appropriate eta for SRF!')
end
