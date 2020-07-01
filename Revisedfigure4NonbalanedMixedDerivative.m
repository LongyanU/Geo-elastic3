clear;
clc

close all


M=7;
% M=4;
AA=zeros(M,M);
b=zeros(M,1);



v=1500;
h=20;
k=linspace(1/50,0.8*pi/h,M);
% k=linspace(1/50,0.43*pi/h,M);

for ii=1:M
    
    for jj=1:M
        AA(ii,jj)=2*cos(jj*k(ii)*h)-2*cos((jj-1)*k(ii)*h);
    end
    b(ii)=-k(ii)^2*h^2;
end


digits(6)
c=AA\b

k=linspace(0,pi/h,100);


for kk=1:5
%  for ii=1:361
%     xita=2*(ii-1)*pi/360;
    xita=(kk-1)*pi/16;
    temp=0;
    for m=1:M
        temp=temp+2*c(m)*sin((m-0.5)*k*h*sin(xita));
    end
    
      temp1=0;
    for m=1:M
        temp1=temp1+2*c(m)*sin((m-0.5)*k*h*cos(xita));
    end
    
     a1=temp.*temp1.*2.*sin((1-0.5)*k*h*cos(xita)).*2.*sin((1-0.5)*k*h*sin(xita))-k.^4*h^4*sin(xita)*cos(xita)*sin(xita)*cos(xita);
%      hold on;plot(a1,'r','LineWidth',2)
     
       if kk==1
        
        figure;plot(a1,'k','LineWidth',2)
    elseif kk==2
        hold on;plot(a1,'m--','LineWidth',2);
    elseif kk==3
        hold on;plot(a1,'r:','LineWidth',2)
    elseif kk==4
        hold on; plot(a1,'b-.','LineWidth',2)
    else
        hold on;plot(a1,'c:.','LineWidth',2)
    end
end


xlabel('percentage of kh')
ylabel('Error')
grid on
legend('¦È=0', '¦È=¦Ð/16','¦È=2¦Ð/16','¦È=3¦Ð/16','¦È=4¦Ð/16')