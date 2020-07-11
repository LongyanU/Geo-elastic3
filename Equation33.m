%Equation 45's correctness

clear;
clc
close all;
nt=1602;
isnap=10;    % snapshot sampling
h=10;
dx=h;
nx=200;
nz=200;
 
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;
 
vp=zeros(nz,nx);
 
vs=vp;
rou=vp;
rou=ones(nz,nx);
f0=70;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^2*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;              % time derivative to obtain gaussian
 
for i=1:200
    for j=1:200
%         vp(i,j)=4908;  %stable
        vp(i,j)=4912;  %unstable
    end
end
 

xs=nz/2;
zs=46;
 
h=dx;
 
%  ******************************************************************
%  txz-------------vx
%  |                |
%  |                |
%  |                |
%  |                |
%  |                |
%  vz---------------txx,tzz
% Virieux 1986 Geophysics
 
tic
 
tic
 
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
 

coeff=[ 1.59906, -0.310692, 0.10345, -0.0398274, 0.0150857, -0.00487876, 0.00104241]; 

1/sqrt(2*sum(abs(coeff)))
 
vs=vp./sqrt(3);
for it=1:nt-2,
    
   Txxx=coeff(1)*( (Txx)-circshift(Txx,[0,1]))+...
        coeff(2)*(circshift(Txx,[0,-1])-circshift(Txx,[0,2]))+...
        coeff(3)*( circshift(Txx,[0,-2])-circshift(Txx,[0,3]))+...
        coeff(4)*( circshift(Txx,[0,-3])-circshift(Txx,[0,4]))+...
        coeff(5)*( circshift(Txx,[0,-4])-circshift(Txx,[0,5]))+...
        coeff(6)*( circshift(Txx,[0,-5])-circshift(Txx,[0,6]))+...
        coeff(7)*( circshift(Txx,[0,-6])-circshift(Txx,[0,7]));
    
    %Txz/z
    %Txzz(i,j)=(Txz(i+1,j)-Txz(i,j));
    Txzz=(circshift(Txz,[ -1])-circshift(Txz,[ 0]));
    
    %Tzz/z
    Tzzz=coeff(1)*( (Tzz)-circshift(Tzz,[1]))+...
        coeff(2)*(circshift(Tzz,[-1])-circshift(Tzz,[2]))+...
        coeff(3)*( circshift(Tzz,[-2])-circshift(Tzz,[3]))+...
        coeff(4)*( circshift(Tzz,[-3])-circshift(Tzz,[4]))+...
        coeff(5)*( circshift(Tzz,[-4])-circshift(Tzz,[5]))+...
        coeff(6)*( circshift(Tzz,[-5])-circshift(Tzz,[6]))+...
        coeff(7)*( circshift(Tzz,[-6])-circshift(Tzz,[7]));
    
    %Txz/x
    %Txzx(i,j)=(Txz(i,j+1)-Txz(i,j));
    Txzx=(circshift(Txz,[0 -1])-circshift(Txz,[0 0]));
    
    Vx=Vx+1./(rou).*dt.*(Txxx+Txzz)/h;
    Vz=Vz+1./(rou).*dt.*(Tzzz+Txzx)/h;
    Vx(zs,xs)=Vx(zs,xs)+src(it);
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    seis_recordVx(it,:)=Vx(zs,:);
    seis_recordVz(it,:)=Vz(zs,:);
    
  Vxz=coeff(1)*(circshift(Vx,[0])-circshift(Vx,[1]))+...
        coeff(2)*(circshift(Vx,[-1])-circshift(Vx,[2]))+...
        coeff(3)*(circshift(Vx,[-2])-circshift(Vx,[3]))+...
        coeff(4)*(circshift(Vx,[-3])-circshift(Vx,[4]))+...
        coeff(5)*(circshift(Vx,[-4])-circshift(Vx,[5]))+...
        coeff(6)*(circshift(Vx,[-5])-circshift(Vx,[6]))+...
        coeff(7)*(circshift(Vx,[-6])-circshift(Vx,[7]));
    
    Vzz=(circshift(Vz,[-1])-circshift(Vz,[0]));
       
        
    Vxx=(circshift(Vx,[0 -1])-circshift(Vx,[0 0]));
      
    
    
    Vzx=coeff(1)*(circshift(Vz,[0 0])-circshift(Vz,[0 1]))+...
        coeff(2)*(circshift(Vz,[0 -1])-circshift(Vz,[0 2]))+...
        coeff(3)*(circshift(Vz,[0 -2])-circshift(Vz,[0 3]))+...
        coeff(4)*(circshift(Vz,[0 -3])-circshift(Vz,[0 4]))+...
        coeff(5)*(circshift(Vz,[0 -4])-circshift(Vz,[0 5]))+...
        coeff(6)*(circshift(Vz,[0 -5])-circshift(Vz,[0 6]))+...
        coeff(7)*(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));
   
    
    Txx=Txx+dt*rou.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*rou.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*rou.*(vs.^2).*(Vxz+Vzx)/h;
    
 
    seis_recordTxz(it,:)=Txz(zs,:);
    
    if rem(it,isnap)== 0,
        imagesc(x,z,Vx), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(Vx))))
        drawnow
    end
 
    
end
 
toc
