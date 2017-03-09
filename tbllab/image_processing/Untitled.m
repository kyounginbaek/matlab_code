close all; clear all; clc;

mu = 1/(4*pi);
gamma = 1/(4*pi);
A=1;
M = 50;
N = 50;

f1 = zeros(M,N);
f2 = zeros(M,N); % 먼저 설정해둬야 속도가 빨라진다

tic % 계산 시간을 측정해준다
for i=1:M
        for j = 1:N
            f(i,j) = A*sin(mu*(i-1) + gamma*(j-1));
        end
end
toc
figure(); surf(f);

tic
r = 0:M-1;
c = 0:N-1;
[C,R] = meshgrid(c,r);
f2 = A*sin(mu*R + gamma*C);
figure(); surf(f2);
toc