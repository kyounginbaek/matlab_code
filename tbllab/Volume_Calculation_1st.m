close all; clear all; clc;

x = 0.7 ; % (cm) ; % 면적의 가로 길이 (for calibration) // N은 1cm, M은 0.7cm
y = 0.9; % (cm); % 면적의 세로 길이 (for calibration) // 랩뷰 라벨 눈금 이미지 - 0.8이면 0.8mm 혹은 0.9면 0.9mm
z = 0.7/(1000-1); % (cm) ; % 사진들 사이 간격의 길이 (for calibraion) // 1000장을 찍었으면 한장한장 사이 길이는 1cm/1000 이나 0.7cm/1000

volume = 0; % 부피를 위한 초기값 지정
  area = 0; % 넓이를 위한 초기값 지정

alpha = 1; % 이미지 시작 넘버
beta = 38; % 받아들인 사진의 개수 (계속 바꿔줘야함!!!)
seta = 49.42857; % hole을 따기 전 사진의 개수

scale = z*(seta/beta); % 사진개수가 다름에 따라 z값이 달라지는걸 보정해주기 위함

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imdata_raw_test = imread('D:\Sample_오리지날\800nm #31 80mJ 00hr_test_01_1.bmp'); % capture 되지 않은 오리지날 이미지를 대입한다
imdata_capture_test = imdata_raw_test(4:408, 29:800, :); % 좌표 추출 (height, width, channel)
figure(); imshow(imdata_raw_test);
figure(); imshow(imdata_capture_test);
[m_test, n_test] = size(imdata_raw_test); % 원본 이미지의 세로(행), 가로(열) 픽셀 수를 받아들인다

for index = alpha : beta
    
    number_pixel = 0; % 픽셀 개수를 세기 위한 초기값 지정
  
    FILEnum = num2str(index);
    FILEname = strcat('D:\Sample_박사님\andy\4th hole(white) 60h 809-852\1 (',FILEnum,').jpg'); % boundary 데이터를 가져온다 -> imfill 필요
    imdata_raw = imread(FILEname); % read image
       
    imdata_bw= im2bw(imdata_raw);
    imdata_fill = imfill(imdata_bw, 'holes'); % imfill을 한다
    [m,n] = size(imdata_fill);
       
    for ii = 1:m
        for jj = 1:n
            if imdata_fill(ii,jj) == 1;
                number_pixel = number_pixel + 1; % 흰면적의 픽셀 개수
               
            end
        end
    end
    
%    number_pixel=length(find(imdata_double(:,:)==1));
    
    if (index == alpha) || (index == beta)
    
        area = (x/n)*(y/m)*number_pixel; % 면적의 넓이 (calibraion)
        volume = volume + (area * scale * 0.5); % 한 이미지당 부피를 계속 더해줌
    
    else 
        
        area = (x/n_test)*(y/m_test)*number_pixel; % 면적의 넓이 (calibraion)
        volume = volume + (area * scale); % 한 이미지당 부피를 계속 더해줌
        
    end
    
    

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

volume % 단위는 cm^3

