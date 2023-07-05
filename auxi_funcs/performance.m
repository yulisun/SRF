function [tp,fp,tn,fn,fplv,fnlv,abfplv,abfnlv,pcc,kappa,imw]=performance(imf,im3)
%����ָ����������
%tp���仯��������fp���������tn��δ�仯��������fn��©����
%��Ҫ��ʽ��tp+fn=Nc; tn+fp=Nu

%im3=double(im3(:,:,3));%�����ᡢ��̫�����ݼ�
im3=double(im3(:,:,1));%�ƺ�һ�š��ƺӶ������ݼ�
imf=double(imf);
[A,B]=size(im3);N=A*B;
Nu=0;Nc=0;
imw=zeros(A,B);%����۲�ͼ
for i=1:A
    for j=1:B
        if im3(i,j)==0
            Nu=Nu+1;
        else
            Nc=Nc+1;
        end
    end
end
im=imf-im3;
fp=0;fn=0;
for i=1:A
    for j=1:B
        if im(i,j)>0
            fp=fp+1;
            imw(i,j)=0;%��ɫ������
        elseif im(i,j)<0
            fn=fn+1;
            imw(i,j)=255;%��ɫ����©��
        else
            imw(i,j)=128;%��ɫ�����޴���
        end
    end
end
imw=uint8(imw);
tp=Nc-fn;tn=Nu-fp;
fplv=fp/N;fnlv=fn/N;
abfplv=fp/Nu;abfnlv=fn/Nc;
pcc=1-fplv-fnlv;

%KAPPAϵ��
pra=(tp+tn)/N;pre=((tp+fp)*(tp+fn)+(fn+tn)*(fp+tn))/(N^2);
kappa=(pra-pre)/(1-pre);