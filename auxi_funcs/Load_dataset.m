if strcmp(dataset,'dataset#1') == 1 % MCD of multispectral VS. multispectral
    load('dataset#1.mat')
    opt.type_t1 = 'optical';
    opt.type_t2 = 'optical';
elseif strcmp(dataset,'dataset#2') == 1 % MCD of multispectral VS. multispectral
    load('dataset#2.mat')
    opt.type_t1 = 'optical';
    opt.type_t2 = 'optical';
elseif strcmp(dataset,'dataset#3') == 1 % MCD of Optical VS. Optical
    load('dataset#3.mat')
    opt.type_t1 = 'optical';
    opt.type_t2 = 'optical';   
elseif strcmp(dataset,'dataset#4') == 1 % MCD of SAR VS. Optical
    load('dataset#4.mat')
    opt.type_t1 = 'sar';
    opt.type_t2 = 'optical';  
elseif strcmp(dataset,'dataset#5') == 1 % MCD of SAR VS. Optical
    load('dataset#5.mat')
    opt.type_t1 = 'sar';
    opt.type_t2 = 'optical';
elseif strcmp(dataset,'dataset#6') == 1 % MCD of Optical VS. SAR
    load('dataset#6.mat')
    opt.type_t1 = 'optical';
    opt.type_t2 = 'sar';
end
%% plot images
figure;
subplot(131);imshow(image_t1);title('imaget1')
subplot(132);imshow(image_t2);title('imaget2')
subplot(133);imshow(Ref_gt,[]);title('Refgt')

