
clear
clc %%%%%%%
close all
nt=4002;    % number of time steps
eps=.6;     % stability
isnap=20;    % snapshot sampling
load('vv')
% 时间已过 1069.186442 秒。
c1=flipud(c);

v=c1;
nx=800;
nx=nx+45*2;
nz=475;
nz=nz+45*2;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=v(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=v(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
        vv(ii+45,jj)=v(ii,800);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end




clear v
v=vv;
vp=v;
vs=vp/sqrt(3);
rho=1000*ones(nz,nx);

dx=10;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;


f0=65;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^2*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian


xs=floor(nx/2);
zs=46;

seis_recordVx=zeros(nt,nx);
seis_recordVz=zeros(nt,nx);
seis_recordTxx=zeros(nt,nx);
seis_recordTzz=zeros(nt,nx);
seis_recordTxz=zeros(nt,nx);

%  ******************************************************************
%  txz-------------vx
%  |                |
%  |                |
%  |                |
%  |                |
%  |                |
%  vz---------------txx,tzz
% Virieux 1986 Geophysics

dx=h;
dz=h;

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


% coeff=[ 1.25219, -0.121591, 0.0331761, -0.0109247, 0.0035117,
% -0.000949881, 0.000177978];    
coeff=[1.22861, -0.102384, 0.0204768, -0.00417893, 0.000689454, -0.0000769225, 0.00000423651]; %Tra FD coefficient

% vp(1:46,:)=vp(1:46,:)*sqrt(2);
% rho(1:46,:)=0.5*10^3;
tic
for it=1:nt-2,
    
    %Txx/x
    Txxx=coeff(1)*( (Txx)-circshift(Txx,[0,1]))+...
        coeff(2)*(circshift(Txx,[0,-1])-circshift(Txx,[0,2]))+...
        coeff(3)*( circshift(Txx,[0,-2])-circshift(Txx,[0,3]))+...
        coeff(4)*( circshift(Txx,[0,-3])-circshift(Txx,[0,4]))+...
        coeff(5)*( circshift(Txx,[0,-4])-circshift(Txx,[0,5]))+...
        coeff(6)*( circshift(Txx,[0,-5])-circshift(Txx,[0,6]))+...
        coeff(7)*( circshift(Txx,[0,-6])-circshift(Txx,[0,7]));
    
    %Txz/z
    %Txzz(i,j)=(Txz(i+1,j)-Txz(i,j));
    Txzz=coeff(1)*(circshift(Txz,[ -1])-circshift(Txz,[ 0]))+...
        coeff(2)*(circshift(Txz,[ -2])-circshift(Txz,[ 1]))+...
        coeff(3)*(circshift(Txz,[ -3])-circshift(Txz,[ 2]))+...
        coeff(4)*(circshift(Txz,[ -4])-circshift(Txz,[ 3]))+...
        coeff(5)*(circshift(Txz,[ -5])-circshift(Txz,[ 4]))+...
        coeff(6)*(circshift(Txz,[ -6])-circshift(Txz,[ 5]))+...
        coeff(7)*(circshift(Txz,[ -7])-circshift(Txz,[ 6]));
        
    
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
    Txzx=coeff(1)*(circshift(Txz,[0 -1])-circshift(Txz,[0 0]))+...
        coeff(2)*(circshift(Txz,[0 -2])-circshift(Txz,[0 1]))+...
        coeff(3)*(circshift(Txz,[0 -3])-circshift(Txz,[0 2]))+...
        coeff(4)*(circshift(Txz,[0 -4])-circshift(Txz,[0 3]))+...
        coeff(5)*(circshift(Txz,[0 -5])-circshift(Txz,[0 4]))+...
        coeff(6)*(circshift(Txz,[0 -6])-circshift(Txz,[0 5]))+...
        coeff(7)*(circshift(Txz,[0 -7])-circshift(Txz,[0 6]));
    
    
    
    Vx=Vx+1./(rho).*dt.*(Txxx+Txzz)/h;
    Vz=Vz+1./(rho).*dt.*(Tzzz+Txzx)/h;
    
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    seismogramVx(it)=Vx(101,400);
    seismogramVz(it)=Vz(101,400);
    
    
    Seismic_Vx(it,:)=Vx(46,:);
    Seismic_Vz(it,:)=Vz(46,:);
    
    
   Vxz=coeff(1)*(circshift(Vx,[0])-circshift(Vx,[1]))+...
        coeff(2)*(circshift(Vx,[-1])-circshift(Vx,[2]))+...
        coeff(3)*(circshift(Vx,[-2])-circshift(Vx,[3]))+...
        coeff(4)*(circshift(Vx,[-3])-circshift(Vx,[4]))+...
        coeff(5)*(circshift(Vx,[-4])-circshift(Vx,[5]))+...
        coeff(6)*(circshift(Vx,[-5])-circshift(Vx,[6]))+...
        coeff(7)*(circshift(Vx,[-6])-circshift(Vx,[7]));
    
    Vzz=coeff(1)*(circshift(Vz,[-1])-circshift(Vz,[0]))+...
        coeff(2)*(circshift(Vz,[-2])-circshift(Vz,[1]))+...
        coeff(3)*(circshift(Vz,[-3])-circshift(Vz,[2]))+...
        coeff(4)*(circshift(Vz,[-4])-circshift(Vz,[3]))+...
        coeff(5)*(circshift(Vz,[-5])-circshift(Vz,[4]))+...
        coeff(6)*(circshift(Vz,[-6])-circshift(Vz,[5]))+...
        coeff(7)*(circshift(Vz,[-7])-circshift(Vz,[6]));
        
    Vxx=coeff(1)*(circshift(Vx,[0 -1])-circshift(Vx,[0 0]))+...
        coeff(2)*(circshift(Vx,[0 -2])-circshift(Vx,[0 1]))+...
        coeff(3)*(circshift(Vx,[0 -3])-circshift(Vx,[0 2]))+...
        coeff(4)*(circshift(Vx,[0 -4])-circshift(Vx,[0 3]))+...
        coeff(5)*(circshift(Vx,[0 -5])-circshift(Vx,[0 4]))+...
        coeff(6)*(circshift(Vx,[0 -6])-circshift(Vx,[0 5]))+...
        coeff(7)*(circshift(Vx,[0 -7])-circshift(Vx,[0 6]));
    
    
    Vzx=coeff(1)*(circshift(Vz,[0 0])-circshift(Vz,[0 1]))+...
        coeff(2)*(circshift(Vz,[0 -1])-circshift(Vz,[0 2]))+...
        coeff(3)*(circshift(Vz,[0 -2])-circshift(Vz,[0 3]))+...
        coeff(4)*(circshift(Vz,[0 -3])-circshift(Vz,[0 4]))+...
        coeff(5)*(circshift(Vz,[0 -4])-circshift(Vz,[0 5]))+...
        coeff(6)*(circshift(Vz,[0 -5])-circshift(Vz,[0 6]))+...
        coeff(7)*(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));
    
    Txx=Txx+dt*rho.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*rho.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*rho.*(vs.^2).*(Vxz+Vzx)/h;
    
    Seismic_Txx(it,:)=Txx(46,:);
    Seismic_Tzz(it,:)=Tzz(46,:);
    Seismic_Txz(it,:)=Txz(46,:);
    
    
    Txx(zs,xs)=Txx(zs,xs)+src(it);
    Tzz(zs,xs)=Tzz(zs,xs)+src(it);
    

    if rem(it,isnap)== 0,
        imagesc(-Txx(46:end,:),[-10^-3 10^-3]), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(Vx))))
        drawnow
    end
end
toc
save('BalancedTra65Hz.mat')