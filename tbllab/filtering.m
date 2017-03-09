a = imread('C:\Users\xnote\Desktop\result_1.jpg');
a = im2double(a);

figure(); imshow(a);

asdf= 1/9;
b = [asdf asdf asdf; asdf asdf asdf; asdf asdf asdf];

[m,n] = size(a);

d = zeros(m+2,n+2);
d(1:m,1:n) = a;

c = zeros(m,n);

for row = 1:m
    for clm = 1:n
        tmp = d(row:row+2,clm:clm+2).*b;
        tmp2 = sum(tmp);
        tmp3 = sum(tmp2);
        c(row,clm) = tmp3;
    end
end

figure(); imshow(c);