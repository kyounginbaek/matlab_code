close all; clear all; clc;

x = 0.7 ; % (cm) ; % ������ ���� ���� (for calibration) // N�� 1cm, M�� 0.7cm
y = 0.9; % (cm); % ������ ���� ���� (for calibration) // ���� �� ���� �̹��� - 0.8�̸� 0.8mm Ȥ�� 0.9�� 0.9mm
z = 0.7/(1000-1); % (cm) ; % ������ ���� ������ ���� (for calibraion) // 1000���� ������� �������� ���� ���̴� 1cm/1000 �̳� 0.7cm/1000

volume = 0; % ���Ǹ� ���� �ʱⰪ ����
  area = 0; % ���̸� ���� �ʱⰪ ����

alpha = 1; % �̹��� ���� �ѹ�
beta = 38; % �޾Ƶ��� ������ ���� (��� �ٲ������!!!)
seta = 49.42857; % hole�� ���� �� ������ ����

scale = z*(seta/beta); % ���������� �ٸ��� ���� z���� �޶����°� �������ֱ� ����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imdata_raw_test = imread('D:\Sample_��������\800nm #31 80mJ 00hr_test_01_1.bmp'); % capture ���� ���� �������� �̹����� �����Ѵ�
imdata_capture_test = imdata_raw_test(4:408, 29:800, :); % ��ǥ ���� (height, width, channel)
figure(); imshow(imdata_raw_test);
figure(); imshow(imdata_capture_test);
[m_test, n_test] = size(imdata_raw_test); % ���� �̹����� ����(��), ����(��) �ȼ� ���� �޾Ƶ��δ�

for index = alpha : beta
    
    number_pixel = 0; % �ȼ� ������ ���� ���� �ʱⰪ ����
  
    FILEnum = num2str(index);
    FILEname = strcat('D:\Sample_�ڻ��\andy\4th hole(white) 60h 809-852\1 (',FILEnum,').jpg'); % boundary �����͸� �����´� -> imfill �ʿ�
    imdata_raw = imread(FILEname); % read image
       
    imdata_bw= im2bw(imdata_raw);
    imdata_fill = imfill(imdata_bw, 'holes'); % imfill�� �Ѵ�
    [m,n] = size(imdata_fill);
       
    for ii = 1:m
        for jj = 1:n
            if imdata_fill(ii,jj) == 1;
                number_pixel = number_pixel + 1; % ������� �ȼ� ����
               
            end
        end
    end
    
%    number_pixel=length(find(imdata_double(:,:)==1));
    
    if (index == alpha) || (index == beta)
    
        area = (x/n)*(y/m)*number_pixel; % ������ ���� (calibraion)
        volume = volume + (area * scale * 0.5); % �� �̹����� ���Ǹ� ��� ������
    
    else 
        
        area = (x/n_test)*(y/m_test)*number_pixel; % ������ ���� (calibraion)
        volume = volume + (area * scale); % �� �̹����� ���Ǹ� ��� ������
        
    end
    
    

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

volume % ������ cm^3

