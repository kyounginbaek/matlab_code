close all; clear all; clc;

x = 7; % ������ ���� ���� (for calibration) // N�� 1cm, M�� 0.7cm
y = 0.9; % ������ ���� ���� (for calibration) // ���� �� ���� �̹��� - 0.8�̸� 0.8mm Ȥ�� 0.9�� 0.9mm%
% z = 10/1000; % ������ ���� ������ ���� (for calibraion) // 1000���� ������� �������� ���� ���̴� 1cm/1000 �̳� 0.7cm/1000

% ������ mm�� �����Ѵ�!

volume = 0; % ���Ǹ� ���� �ʱⰪ ����
  area = 0; % ���̸� ���� �ʱⰪ ����
  
alpha = 1; % �̹��� ���� �ѹ�
beta = 8; % �޾Ƶ��� ������ ���� (��� �ٲ������!!!)
scale = 0.1; % volume�� ���ϱ� ���� plot�� �Ҷ� x���� �� �� ����ȭ �ϱ����� scale

A = zeros(1, beta-alpha+1); % area plot�� ���� ��� A ����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index = alpha : beta
    
    number_pixel = 0; % �ȼ� ������ ���� ���� �ʱⰪ ����
  
    
    FILEnum = num2str(index);
    FILEname = strcat('C:\Users\xnote\Desktop\TBL LAB\Sample - label max\1 (',FILEnum,').bmp');
    imdata_raw = imread(FILEname); % read image
    [m,n] = size(imdata_raw); % ����, ���� �ȼ� ���� �޾Ƶ��δ�
    
    imdata_double = im2double(imdata_raw);
    
    for ii = 1:m 
        for jj = 1:n
            if imdata_double(ii,jj) == 1;
                number_pixel = number_pixel + 1; % ������� �ȼ� ����
               
            end
        end
    end
    
    area = (x/n)*(y/m)*number_pixel; % ������ ���� (calibraion)
    
    A(index) = area;

end

    x = alpha : beta; % plot�� ���� x���� ���� ����
    y = A(x); % plot�� ���� y�� ����
    xx = alpha : scale : beta; % plot�� x���� ������ �ɰ���
    yy = spline(x,y,xx); 
    plot(x,y,'o',xx,yy)
    
    
for index = alpha : 1 : (beta-1)*(1/scale)+1
     
    volume = volume + yy(index); % �ɰ��� area���� ��� ���ؼ� volume�� ���Ѵ�
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
