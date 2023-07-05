function [regression_t1,regression_t2,DI_t1,DI_t2,CM] = SRF_main(image_t1,image_t2,opt)
%-------------  Preprocessing---------------%
if strcmp(opt.type_t1 ,'sar') 
    image_t1 = image_normlized(image_t1,opt.type_t1);
end
if strcmp(opt.type_t2,'sar') 
    image_t2 = image_normlized(image_t2,opt.type_t2);
end
image_t1 = double(image_t1);
image_t2 = double(image_t2);
h = fspecial('average',5);
image_t1 = imfilter(image_t1, h,'symmetric');
image_t2 = imfilter(image_t2, h,'symmetric');
[Cosup,~] = GMMSP_Cosegmentation(image_t1,image_t2,opt.Ns);
[t1_feature,t2_feature,norm_par] = MMfeature_extraction(Cosup,image_t1,image_t2); 

%-------------  Structure representation---------------%
opt.kmax =round(sqrt(size(t1_feature,2))*1);
[ADJ_t1] = AdaptiveStructureGraph(t1_feature',t1_feature',opt.kmax);
[ADJ_t2] = AdaptiveStructureGraph(t2_feature',t2_feature',opt.kmax);
ADJ_t1 = ADJ_t1-diag(diag(ADJ_t1));
ADJ_t2 = ADJ_t2-diag(diag(ADJ_t2));
ADJ_fuse = min (ADJ_t1,ADJ_t2);
Lx = LaplacianMatrix(ADJ_t1);
Ly = LaplacianMatrix(ADJ_t2);
Lf = LaplacianMatrix(ADJ_fuse);

%-------------  Structural regression fusion---------------%
[Zx, delt_t1,Zy, delt_t2,~] = Structural_regression_fusion(t1_feature,t2_feature,Lx,Ly,Lf,opt);
[regression_t1,~,~] = suplabel2ImFeature(Cosup,Zy,size(image_t1,3));% t2--->t1
regression_t1 = DenormImage(regression_t1,norm_par(1:size(image_t1,3)));
DI_t1  = suplabel2DI(Cosup,sum(delt_t2.^2,1));

[regression_t2,~,~] = suplabel2ImFeature(Cosup,Zx,size(image_t2,3));% t1--->t2
regression_t2 = DenormImage(regression_t2,norm_par(size(image_t1,3)+1:end));
DI_t2  = suplabel2DI(Cosup,sum(delt_t1.^2,1));
%-------------  Change extraction---------------%
fx = sqrt(sum(delt_t1.^2,1));
fy = sqrt(sum(delt_t2.^2,1));
[CM,~] = MRF_CoSegmentation(Cosup,opt.alfa,t1_feature,t2_feature,fx,fy);
%%


