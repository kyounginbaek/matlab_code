close all; clear all; clc;

x = 7; % 면적의 가로 길이 (for calibration) // N은 1cm, M은 0.7cm
y = 0.9; % 면적의 세로 길이 (for calibration) // 랩뷰 라벨 눈금 이미지 - 0.8이면 0.8mm 혹은 0.9면 0.9mm%
% z = 10/1000; % 사진들 사이 간격의 길이 (for calibraion) // 1000장을 찍었으면 한장한장 사이 길이는 1cm/1000 이나 0.7cm/1000

% 단위는 mm로 통일한다!

volume = 0; % 부피를 위한 초기값 지정
  area = 0; % 넓이를 위한 초기값 지정
  
alpha = 1; % 이미지 시작 넘버
beta = 8; % 받아들인 사진의 개수 (계속 바꿔줘야함!!!)
scale = 0.1; % volume을 구하기 위한 plot을 할때 x축을 좀 더 세분화 하기위한 scale

A = zeros(1, beta-alpha+1); % area plot을 위한 행렬 A 생성

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index = alpha : beta
    
    number_pixel = 0; % 픽셀 개수를 세기 위한 초기값 지정
  
    
    FILEnum = num2str(index);
    FILEname = strcat('C:\Users\xnote\Desktop\TBL LAB\Sample - label max\1 (',FILEnum,').bmp');
    imdata_raw = imread(FILEname); % read image
    [m,n] = size(imdata_raw); % 세로, 가로 픽셀 수를 받아들인다
    
    imdata_double = im2double(imdata_raw);
    
    for ii = 1:m 
        for jj = 1:n
            if imdata_double(ii,jj) == 1;
                number_pixel = number_pixel + 1; % 흰면적의 픽셀 개수
               
            end
        end
    end
    
    area = (x/n)*(y/m)*number_pixel; % 면적의 넓이 (calibraion)
    
    A(index) = area;

end

    x = alpha : beta; % plot을 위한 x값의 범위 설정
    y = A(x); % plot을 위한 y값 설정
    xx = alpha : scale : beta; % plot의 x값의 범위를 쪼갠다
    yy = spline(x,y,xx); 
    plot(x,y,'o',xx,yy)
    
    
for index = alpha : 1 : (beta-1)*(1/scale)+1
     
    volume = volume + yy(index); % 쪼개진 area들을 계속 더해서 volume을 구한다
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
