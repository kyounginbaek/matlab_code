close all; clear all; clc;

alpha = 1; % �̹��� ���� �ѹ�
beta = 8; % �޾Ƶ��� ������ ���� (��� �ٲ������!!!)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index1 = alpha : beta;
    
    number_pixel = 0; % �ȼ� ������ ���� ���� �ʱⰪ ����
    
    FILEnum = num2str(index1);
    FILEname = strcat('D:\TLB Lab (���� - MATLAB)\���� Sample (Volume Calculation)\1 (',FILEnum,').bmp');
    imdata_raw = imread(FILEname); % read image
    [m,n] = size(imdata_raw); % ����,���� �ȼ� ���� �޾Ƶ��δ�
    
    imdata_double = im2double (imdata_raw);
    
    for ii = 1:m
        for jj = 1:n
            if imdata_double(ii,jj) < 1;
                imdata_double(ii,jj) = 0;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [l,n]=bwlabel(imdata_double);

    max_points=find(l==1);

    for index2 = 2:n
        
        points = find(l==index2);
        
        if(length(points) > length(max_points))
            max_points = points;
            
        end
        
    end
    
    
    [r,c]=size(imdata_double);
    label_result=zeros(r,c);
    
    for index3 = 1:length(max_points)
        
        label_result(max_points(index3)) = 1;
        
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = strcat('result',num2str(index1),'.bmp');
imwrite(label_result,name);

close all;

end
