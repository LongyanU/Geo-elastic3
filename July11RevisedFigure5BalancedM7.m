% Please first run Garvin2(0.5,0.55) to get GarvinResult.mat
% and then run this program so as to compare with analytic results
% 时间已过 915.457025 秒。
clear;
clc
close all;

% nt=6346;
nt=3346*1;
nt=5502;
isnap=500;    % snapshot sampling

dx=10;
h=10;
%nx=430;
nx=900;
%nz=220;
nz=400;

x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
% dt=0.003; % calculate time step from stability criterion
dt=0.001; % calculate time step from stability criterion
tau=dt;


t=(1:nt)*dt;


f0=14;  %f0=1/tp %the peak frequency of the seismic source
t0=1.5/f0;     %ts

% src=10^7*exp(-(pi*f0*(t-t0)).^2);              % source time function
% src=-diff((src))/dx^2;				% time derivative to obtain gaussian
src=10^7*((pi*f0*(t-t0)).^2-0.5).*exp(-(pi*f0*(t-t0)).^2);



for i=1:nz
    for j=1:nx
        %vp(i,j)=4000;
        vp(i,j)=1732.10;
    end
end


for i=1:nz
    for j=1:nx
        %vs(i,j)=2310;
        vs(i,j)=1000.0;
    end
end

for i=1:nz
    for j=1:nx
        %rho(i,j)=2.5*10^6;
        rho(i,j)=1.0*10^3;
    end
end

vp(1:46*1,:)=1000.0*1^0.5; 
rho(1:46*1,:)=0.5*10^3; 

p=zeros([nz nx]); Vx=p; Vz=p;
Txxx=p;
Txzz=p;
Tzzz=p;
Txzx=p;

Txx=p;
Txz=p;
Tzz=p;


Vxx=p;
Vzz=p;
Vxz=p;
Vzx=p;

coeff=[ 1.24558, -0.115879, 0.0289305, -0.00826984, 0.00215709, -0.000433658, 0.0000500198];
ixs=200*1;
% izs= 46*1+1+200; %%%%%%
izs= 46+1; %%%%%%

Seismic_Vx=zeros(nt,nx);
Seismic_Vz=zeros(nt,nx);
Seismic_Txx=zeros(nt,nx);
Seismic_Tzz=zeros(nt,nx);
Seismic_Txz=zeros(nt,nx);

tic
for it=1:nt-2,
    Txxx=coeff(1)*( (Txx)-circshift(Txx,[0,-1]))+...
        coeff(2)*(circshift(Txx,[0,1])-circshift(Txx,[0,-2]))+...
        coeff(3)*( circshift(Txx,[0,2])-circshift(Txx,[0,-3]))+...
        coeff(4)*( circshift(Txx,[0,3])-circshift(Txx,[0,-4]))+...
        coeff(5)*( circshift(Txx,[0,4])-circshift(Txx,[0,-5]))+...
        coeff(6)*( circshift(Txx,[0,5])-circshift(Txx,[0,-6]))+...
        coeff(7)*( circshift(Txx,[0,6])-circshift(Txx,[0,-7]));
    
    %Txz/z
    %Txzz(i,j)=(Txz(i+1,j)-Txz(i,j));
    Txzz=coeff(1)*(circshift(Txz,[ 1])-circshift(Txz,[ 0]))+...
        coeff(2)*(circshift(Txz,[ 2])-circshift(Txz,[ -1]))+...
        coeff(3)*(circshift(Txz,[ 3])-circshift(Txz,[ -2]))+...
        coeff(4)*(circshift(Txz,[ 4])-circshift(Txz,[ -3]))+...
        coeff(5)*(circshift(Txz,[ 5])-circshift(Txz,[ -4]))+...
        coeff(6)*(circshift(Txz,[ 6])-circshift(Txz,[ -5]))+...
        coeff(7)*(circshift(Txz,[ 7])-circshift(Txz,[ -6]));
    
    %Tzz/z
    Tzzz=coeff(1)*( (Tzz)-circshift(Tzz,[-1]))+...
        coeff(2)*(circshift(Tzz,[1])-circshift(Tzz,[-2]))+...
        coeff(3)*( circshift(Tzz,[2])-circshift(Tzz,[-3]))+...
        coeff(4)*( circshift(Tzz,[3])-circshift(Tzz,[-4]))+...
        coeff(5)*( circshift(Tzz,[4])-circshift(Tzz,[-5]))+...
        coeff(6)*( circshift(Tzz,[5])-circshift(Tzz,[-6]))+...
        coeff(7)*( circshift(Tzz,[6])-circshift(Tzz,[-7]));
    
    %Txz/x
    %Txzx(i,j)=(Txz(i,j+1)-Txz(i,j));
    Txzx=coeff(1)*(circshift(Txz,[0 1])-Txz)+...
        coeff(2)*(circshift(Txz,[0 2])-circshift(Txz,[0 -1]))+...
        coeff(3)*(circshift(Txz,[0 3])-circshift(Txz,[0 -2]))+...
        coeff(4)*(circshift(Txz,[0 4])-circshift(Txz,[0 -3]))+...
        coeff(5)*(circshift(Txz,[0 5])-circshift(Txz,[0 -4]))+...
        coeff(6)*(circshift(Txz,[0 6])-circshift(Txz,[0 -5]))+...
        coeff(7)*(circshift(Txz,[0 7])-circshift(Txz,[0 -6]));
    
    
    Vx=Vx+1./(rho).*dt.*(Txxx+Txzz)/h;
    Vz=Vz+1./(rho).*dt.*(Tzzz+Txzx)/h;
    
    Vx(izs,ixs)=Vx(izs,ixs)+src(it);
    Vz(izs,ixs)=Vz(izs,ixs)+src(it);
    
%     [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,85,85,0.005);
         [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
%     if it==500
%         snapshotVx1=Vx;
%         snapshotVz1=Vz;
%     elseif it==1000
%         snapshotVx2=Vx;
%         snapshotVz2=Vz;
%     elseif it==1500
%         snapshotVx3=Vx;
%         snapshotVz3=Vz;
%     elseif it==2000
%         snapshotVx4=Vx;
%         snapshotVz4=Vz;
%     elseif it==2500
%         snapshotVx5=Vx;
%         snapshotVz5=Vz;
%     elseif it==3000
%         snapshotVx6=Vx;
%         snapshotVz6=Vz;
%     end
%     
    seismogramVx(it)=Vx(101*1+1,250*1);
    seismogramVz(it)=Vz(101*1+1,250*1);
    
    
    seismogramVx2(it)=Vx(101*1+1,300*1);
    seismogramVz2(it)=Vz(101*1+1,300*1);
    
    seismogramVx3(it)=Vx(101*1+1,400*1);
    seismogramVz3(it)=Vz(101*1+1,400*1);
    
    seismogramVx4(it)=Vx(101*1+1,500*1);
    seismogramVz4(it)=Vz(101*1+1,500*1);
    
    seismogramVx5(it)=Vx(101*1+1,600*1);
    seismogramVz5(it)=Vz(101*1+1,600*1);
    
    seismogramVx6(it)=Vx(101*1+1,650*1);
    seismogramVz6(it)=Vz(101*1+1,650*1);
    
    seismogramVx7(it)=Vx(150*1+1,400*1);
    seismogramVz7(it)=Vz(150*1+1,400*1);
    
    seismogramV8(it)=Vx(200*1+1,400*1);
    seismogramVz8(it)=Vz(200*1+1,400*1);
    
    seismogramVx9(it)=Vx(250*1+1,400*1);
    seismogramVz9(it)=Vz(250*1+1,400*1);
    
    %%%
    seismogramVx10(it)=Vx(101*1+1,700*1);
    seismogramVz10(it)=Vz(101*1+1,700*1);
    
    seismogramVx11(it)=Vx(101*1+1,750*1);
    seismogramVz11(it)=Vz(101*1+1,750*1);
    
    
    seismogramVx12(it)=Vx(150*1+1,750*1);
    seismogramVz12(it)=Vz(150*1+1,750*1);
    
    seismogramVx13(it)=Vx(200*1+1,750*1);
    seismogramVz13(it)=Vz(200*1+1,750*1);
    
    seismogramVx14(it)=Vx(250*1+1,750*1);
    seismogramVz14(it)=Vz(250*1+1,750*1);
    
    seismogramVx15(it)=Vx(300*1+1,750*1);
    seismogramVz15(it)=Vz(300*1+1,750*1);
    
    %%%%%
    

    
Vxz=coeff(1)*(circshift(Vx,[0])-circshift(Vx,[-1]))+...
        coeff(2)*(circshift(Vx,[1])-circshift(Vx,[-2]))+...
        coeff(3)*(circshift(Vx,[2])-circshift(Vx,[-3]))+...
        coeff(4)*(circshift(Vx,[3])-circshift(Vx,[-4]))+...
        coeff(5)*(circshift(Vx,[4])-circshift(Vx,[-5]))+...
        coeff(6)*(circshift(Vx,[5])-circshift(Vx,[-6]))+...
        coeff(7)*(circshift(Vx,[6])-circshift(Vx,[-7]));
    
    Vzz=coeff(1)*(circshift(Vz,[1])-circshift(Vz,[0]))+...
        coeff(2)*(circshift(Vz,[2])-circshift(Vz,[-1]))+...
        coeff(3)*(circshift(Vz,[3])-circshift(Vz,[-2]))+...
        coeff(4)*(circshift(Vz,[4])-circshift(Vz,[-3]))+...
        coeff(5)*(circshift(Vz,[5])-circshift(Vz,[-4]))+...
        coeff(6)*(circshift(Vz,[6])-circshift(Vz,[-5]))+...
        coeff(7)*(circshift(Vz,[7])-circshift(Vz,[-6]));
        
    Vxx=coeff(1)*(circshift(Vx,[0 1])-circshift(Vx,[0 0]))+...
        coeff(2)*(circshift(Vx,[0 2])-circshift(Vx,[0 -1]))+...
        coeff(3)*(circshift(Vx,[0 3])-circshift(Vx,[0 -2]))+...
        coeff(4)*(circshift(Vx,[0 4])-circshift(Vx,[0 -3]))+...
        coeff(5)*(circshift(Vx,[0 5])-circshift(Vx,[0 -4]))+...
        coeff(6)*(circshift(Vx,[0 6])-circshift(Vx,[0 -5]))+...
        coeff(7)*(circshift(Vx,[0 7])-circshift(Vx,[0 -6]));

        
    
    Vzx=coeff(1)*(circshift(Vz,[0 0])-circshift(Vz,[0 -1]))+...
        coeff(2)*(circshift(Vz,[0 1])-circshift(Vz,[0 -2]))+...
        coeff(3)*(circshift(Vz,[0 2])-circshift(Vz,[0 -3]))+...
        coeff(4)*(circshift(Vz,[0 3])-circshift(Vz,[0 -4]))+...
        coeff(5)*(circshift(Vz,[0 4])-circshift(Vz,[0 -5]))+...
        coeff(6)*(circshift(Vz,[0 5])-circshift(Vz,[0 -6]))+...
        coeff(7)*(circshift(Vz,[0 6])-circshift(Vz,[0 -7]));
    
    %         Txx(47:end,:)=Txx(47:end,:)+dt/h*rho(47:end,:).*(vp(47:end,:).^2.*Vxx(47:end,:)+(vp(47:end,:).^2-2*vs(47:end,:).^2).*Vzz(47:end,:));
    Txx=Txx+dt*rho.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*rho.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*rho.*(vs.^2).*(Vxz+Vzx)/h;
    
    %         Txx(1:46,:)=Txx(1:46,:)+dt*4*rho(1:46,:).*vs(1:46,:).^2./vp(1:46,:).^2.*((vp(1:46,:).^2-vs(1:46,:).^2).*Vxx(1:46,:))/h;
    
    

    
    
%     Txx(izs,ixs)=Txx(izs,ixs)+src(it);
%     Tzz(izs,ixs)=Tzz(izs,ixs)+src(it);
    
    
    Tzz(1:46*1,:)=0;
%     Txz(1:46*1,:)=0;%%

    if rem(it,isnap)== 0,
        imagesc(-Vx(46*1:end,:),[-10^5 10^5]), axis equal
        hold on; plot(200*1,2,'r*','linewidth',2)
        
        hold on; plot(250*1,55*1,'bo','linewidth',2)
        hold on; plot(300*1,55*1,'bo','linewidth',2)
        hold on; plot(400*1,55*1,'bo','linewidth',2)
        hold on; plot(500*1,55*1,'bo','linewidth',2)
        hold on; plot(600*1,55*1,'bo','linewidth',2)
        hold on; plot(650*1,55*1,'bo','linewidth',2)
        
        hold on; plot(400*1,(150-46)*1,'bo','linewidth',2)
        hold on; plot(400*1,(200-46)*1,'bo','linewidth',2)
        hold on; plot(400*1,(250-46)*1,'bo','linewidth',2)
        
        hold on; plot(700*1,55*1,'bo','linewidth',2)
        hold on; plot(750*1,55*1,'bo','linewidth',2)
        hold on; plot(750*1,(150-46)*1,'bo','linewidth',2)
        
        hold on; plot(750*1,(200-46)*1,'bo','linewidth',2)
        hold on; plot(750*1,(250-46)*1,'bo','linewidth',2)
        hold on; plot(750*1,(300-46)*1,'bo','linewidth',2)
        
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(Vx))))
        drawnow
    end
end
toc


% save('NonbalancedM7July10version2.mat')
save('balancedM7.mat')
% save('NonbalancedSmallGrid.mat')
figure;imagesc(vs(46*1:end,:))
% hold on; plot(200*1,200,'r*','linewidth',2)
hold on; plot(200*1,2,'r*','linewidth',2)
% izs= 46*1+1+200; %%%%%%
% ixs=200*1;
hold on; plot(250*1,55*1,'bo','linewidth',2)
hold on; plot(300*1,55*1,'bo','linewidth',2)
hold on; plot(400*1,55*1,'bo','linewidth',2)
hold on; plot(500*1,55*1,'bo','linewidth',2)
hold on; plot(600*1,55*1,'bo','linewidth',2)
hold on; plot(650*1,55*1,'bo','linewidth',2)

hold on; plot(400*1,(150-46)*1,'bo','linewidth',2)
hold on; plot(400*1,(200-46)*1,'bo','linewidth',2)
hold on; plot(400*1,(250-46)*1,'bo','linewidth',2)

hold on; plot(700*1,55*1,'bo','linewidth',2)
hold on; plot(750*1,55*1,'bo','linewidth',2)
hold on; plot(750*1,(150-46)*1,'bo','linewidth',2)

hold on; plot(750*1,(200-46)*1,'bo','linewidth',2)
hold on; plot(750*1,(250-46)*1,'bo','linewidth',2)
hold on; plot(750*1,(300-46)*1,'bo','linewidth',2)

grid on
legend('the red is the seismic location','The blue is the receiver location' );
xlabel('x/dx')
ylabel('z/dz')
