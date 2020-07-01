clear;
clc

% close all
% [ 1.46951, -0.117903] 已经取得了很不错的结果

% 程序解释 p(1)+p(-1)-2p(0) =2cos(kh)-2   =   A
% 程序解释 p(2)+p(-2)-2p(0) =2cos(2kh)-2  =B
% coeff(1)*A + coeff(2)*B= -k^2*h^2;
% coeff是需要求解的系数


M=7;
AA=zeros(M,M);
b=zeros(M,1);


v=1500;
h=20;
k=linspace(1/50,0.61*pi/h,M);

for ii=1:M
    
    for jj=1:M
        AA(ii,jj)=2*sin((jj-0.5)*k(ii)*h);
    end
    b(ii)=k(ii)*h;
end


digits(6)
c=AA\b

k=linspace(0,pi/h,100);

temp=0;
for m=1:M
    temp=temp+2*c(m)*sin((m-0.5)*k*h);
end

a1=temp.^2-k.^2*h^2;
figure;plot(a1,'k','LineWidth',2)



xlabel('percentage of kh')
ylabel('error')
grid on
axis([0 100 -5*10^-4  2*10^-4])