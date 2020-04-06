% first run Figure11aStaggeredFDTra.m,figure11bStaggeredFD.m,figure11cKspace to get seismic records,
% then run figure11compareSeimogramsVz.m, Figure11CompareSeis_recordVz.m; figure12compareSeimogramsTxx.m,figure12CompareSeis_recordTxx
% to get the figures in figure 11 and figure 12.
% This is only for the convenient of the reviewers.

% 时间已过 2355.776368 秒。  Dec 19
% 时间已过 3780.086018 秒。
% 时间已过 4337.352889 秒。
clear
clc %%%%%%%
close all
nt=4002;    % number of time steps
eps=.6;     % stability
isnap=20;    % snapshot sampling
load('vv')

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
% coeff=[1.22861, -0.102384, 0.0204768, -0.00417893, 0.000689454, -0.0000769225, 0.00000423651]; %Tra FD coefficient

coeff=[ 1.25885, -0.12771, 0.0383129, -0.014856, 0.00622176, -0.0026224, 0.0010729, -0.000416827, 0.000151442, -0.0000508417,...
    0.000015607, -0.00000433692, 0.00000107964, -2.38015e-7, 4.58484e-8, -7.5909e-9, 1.05768e-9, -1.20546e-10, ...
    1.07866e-11, -7.10391e-13, 3.06084e-14, -6.47146e-16];

coeff=[ 1.26267, -0.131246, 0.0413424, -0.017258, 0.00798353, -0.00381739, 0.00182211, -0.000850754, 0.000383467, -0.0001653, 0.0000676569, -0.0000261375, 0.0000094812, -0.00000321363, 0.00000101297, -2.95492e-7, 7.93619e-8, -1.95141e-8, 4.36538e-9, -8.82049e-10, 1.59619e-10, -2.56087e-11, 3.59738e-12, -4.3555e-13, 4.45245e-14, -3.73644e-15, 2.47126e-16, -1.20779e-17, 3.87765e-19, -6.13426e-21]
 


% vp(1:46,:)=vp(1:46,:)*sqrt(2);
% rho(1:46,:)=0.5*10^3;
tic
for it=1:nt-2,
    
    %Txx/x
    %     Txxx=coeff(1)*( (Txx)-circshift(Txx,[0,1]))+...
    %         coeff(2)*(circshift(Txx,[0,-1])-circshift(Txx,[0,2]))+...
    %         coeff(3)*( circshift(Txx,[0,-2])-circshift(Txx,[0,3]))+...
    %         coeff(4)*( circshift(Txx,[0,-3])-circshift(Txx,[0,4]))+...
    %         coeff(5)*( circshift(Txx,[0,-4])-circshift(Txx,[0,5]))+...
    %         coeff(6)*( circshift(Txx,[0,-5])-circshift(Txx,[0,6]))+...
    %         coeff(7)*( circshift(Txx,[0,-6])-circshift(Txx,[0,7]))+...
    %         coeff(8)*( circshift(Txx,[0,-7])-circshift(Txx,[0,8]))+...
    %         coeff(9)*( circshift(Txx,[0,-8])-circshift(Txx,[0,9]))+...
    %         coeff(10)*( circshift(Txx,[0,-9])-circshift(Txx,[0,10]))+...
    %         coeff(11)*( circshift(Txx,[0,-10])-circshift(Txx,[0,11]))+...
    %         coeff(12)*( circshift(Txx,[0,-11])-circshift(Txx,[0,12]))+...
    %         coeff(13)*( circshift(Txx,[0,-12])-circshift(Txx,[0,13]))+...
    %         coeff(14)*( circshift(Txx,[0,-13])-circshift(Txx,[0,14]))+...
    %         coeff(15)*( circshift(Txx,[0,-14])-circshift(Txx,[0,15]))+...
    %         coeff(16)*( circshift(Txx,[0,-15])-circshift(Txx,[0,16]))+...
    %         coeff(17)*( circshift(Txx,[0,-16])-circshift(Txx,[0,17]))+...
    %         coeff(18)*( circshift(Txx,[0,-17])-circshift(Txx,[0,18]))+...
    %         coeff(19)*( circshift(Txx,[0,-18])-circshift(Txx,[0,19]))+...
    %         coeff(20)*( circshift(Txx,[0,-19])-circshift(Txx,[0,20]))+...
    %         coeff(21)*( circshift(Txx,[0,-20])-circshift(Txx,[0,21]))+...
    %         coeff(22)*( circshift(Txx,[0,-21])-circshift(Txx,[0,22]));
    
    
    Txxx=coeff(1)*( (Txx)-circshift(Txx,[0,1]));
    for iii=2:length(coeff)
        Txxx=Txxx+coeff(iii)*(circshift(Txx,[0,-iii+1])-circshift(Txx,[0,iii]));
    end
    
    
    %Txz/z
    %Txzz(i,j)=(Txz(i+1,j)-Txz(i,j));
    Txzz=coeff(1)*(circshift(Txz,[ -1])-circshift(Txz,[ 0]));
    for iii=2:length(coeff)
        Txzz=Txzz+coeff(iii)*(circshift(Txz,[ -iii])-circshift(Txz,[ iii-1]));
    end
    
    %     Txzz=coeff(1)*(circshift(Txz,[ -1])-circshift(Txz,[ 0]))+...
    %         coeff(2)*(circshift(Txz,[ -2])-circshift(Txz,[ 1]))+...
    %         coeff(3)*(circshift(Txz,[ -3])-circshift(Txz,[ 2]))+...
    %         coeff(4)*(circshift(Txz,[ -4])-circshift(Txz,[ 3]))+...
    %         coeff(5)*(circshift(Txz,[ -5])-circshift(Txz,[ 4]))+...
    %         coeff(6)*(circshift(Txz,[ -6])-circshift(Txz,[ 5]))+...
    %         coeff(7)*(circshift(Txz,[ -7])-circshift(Txz,[ 6]));
    
    
    
    
    %Tzz/z
    
    Tzzz=coeff(1)*( (Tzz)-circshift(Tzz,[1]));
    for iii=2:length(coeff)
        Tzzz=Tzzz+coeff(iii)*(circshift(Tzz,[-iii+1])-circshift(Tzz,[iii]));
    end
    %
    %     Tzzz=coeff(1)*( (Tzz)-circshift(Tzz,[1]))+...
    %         coeff(2)*(circshift(Tzz,[-1])-circshift(Tzz,[2]))+...
    %         coeff(3)*( circshift(Tzz,[-2])-circshift(Tzz,[3]))+...
    %         coeff(4)*( circshift(Tzz,[-3])-circshift(Tzz,[4]))+...
    %         coeff(5)*( circshift(Tzz,[-4])-circshift(Tzz,[5]))+...
    %         coeff(6)*( circshift(Tzz,[-5])-circshift(Tzz,[6]))+...
    %         coeff(7)*( circshift(Tzz,[-6])-circshift(Tzz,[7]));
    
    %Txz/x
    %Txzx(i,j)=(Txz(i,j+1)-Txz(i,j));
    
    Txzx=coeff(1)*(circshift(Txz,[0 -1])-circshift(Txz,[0 0]));
    for iii=2:length(coeff)
        Txzx=Txzx+ coeff(iii)*(circshift(Txz,[0 -iii])-circshift(Txz,[0 iii-1]));
    end
    %     Txzx=coeff(1)*(circshift(Txz,[0 -1])-circshift(Txz,[0 0]))+...
    %         coeff(2)*(circshift(Txz,[0 -2])-circshift(Txz,[0 1]))+...
    %         coeff(3)*(circshift(Txz,[0 -3])-circshift(Txz,[0 2]))+...
    %         coeff(4)*(circshift(Txz,[0 -4])-circshift(Txz,[0 3]))+...
    %         coeff(5)*(circshift(Txz,[0 -5])-circshift(Txz,[0 4]))+...
    %         coeff(6)*(circshift(Txz,[0 -6])-circshift(Txz,[0 5]))+...
    %         coeff(7)*(circshift(Txz,[0 -7])-circshift(Txz,[0 6]));
    
    
    
    Vx=Vx+1./(rho).*dt.*(Txxx+Txzz)/h;
    Vz=Vz+1./(rho).*dt.*(Tzzz+Txzx)/h;
    
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    seismogramVx(it)=Vx(101,400);
    seismogramVz(it)=Vz(101,400);
    
    
    
    Seismic_Vx(it,:)=Vx(46,:);
    Seismic_Vz(it,:)=Vz(46,:);
    
    
    %     Vxz=coeff(1)*(circshift(Vx,[0])-circshift(Vx,[1]))+...
    %         coeff(2)*(circshift(Vx,[-1])-circshift(Vx,[2]))+...
    %         coeff(3)*(circshift(Vx,[-2])-circshift(Vx,[3]))+...
    %         coeff(4)*(circshift(Vx,[-3])-circshift(Vx,[4]))+...
    %         coeff(5)*(circshift(Vx,[-4])-circshift(Vx,[5]))+...
    %         coeff(6)*(circshift(Vx,[-5])-circshift(Vx,[6]))+...
    %         coeff(7)*(circshift(Vx,[-6])-circshift(Vx,[7]));
    Vxz=coeff(1)*(circshift(Vx,[0])-circshift(Vx,[1]));
    for iii=2:length(coeff)
        Vxz=Vxz+coeff(iii)*(circshift(Vx,[-iii+1])-circshift(Vx,[iii]));
    end
    
    
    %     Vzz=coeff(1)*(circshift(Vz,[-1])-circshift(Vz,[0]))+...
    %         coeff(2)*(circshift(Vz,[-2])-circshift(Vz,[1]))+...
    %         coeff(3)*(circshift(Vz,[-3])-circshift(Vz,[2]))+...
    %         coeff(4)*(circshift(Vz,[-4])-circshift(Vz,[3]))+...
    %         coeff(5)*(circshift(Vz,[-5])-circshift(Vz,[4]))+...
    %         coeff(6)*(circshift(Vz,[-6])-circshift(Vz,[5]))+...
    %         coeff(7)*(circshift(Vz,[-7])-circshift(Vz,[6]));
    Vzz=coeff(1)*(circshift(Vz,[-1])-circshift(Vz,[0]));
    for iii=2:length(coeff)
        Vzz=Vzz+ coeff(iii)*(circshift(Vz,[-iii])-circshift(Vz,[iii-1]));
    end
    
    
    %     Vxx=coeff(1)*(circshift(Vx,[0 -1])-circshift(Vx,[0 0]))+...
    %         coeff(2)*(circshift(Vx,[0 -2])-circshift(Vx,[0 1]))+...
    %         coeff(3)*(circshift(Vx,[0 -3])-circshift(Vx,[0 2]))+...
    %         coeff(4)*(circshift(Vx,[0 -4])-circshift(Vx,[0 3]))+...
    %         coeff(5)*(circshift(Vx,[0 -5])-circshift(Vx,[0 4]))+...
    %         coeff(6)*(circshift(Vx,[0 -6])-circshift(Vx,[0 5]))+...
    %         coeff(7)*(circshift(Vx,[0 -7])-circshift(Vx,[0 6]));
    Vxx=coeff(1)*(circshift(Vx,[0 -1])-circshift(Vx,[0 0]));
    for iii=2:length(coeff)
        Vxx=Vxx+ coeff(iii)*(circshift(Vx,[0 -iii])-circshift(Vx,[0 iii-1]));
    end
    
    %     Vzx=coeff(1)*(circshift(Vz,[0 0])-circshift(Vz,[0 1]))+...
    %         coeff(2)*(circshift(Vz,[0 -1])-circshift(Vz,[0 2]))+...
    %         coeff(3)*(circshift(Vz,[0 -2])-circshift(Vz,[0 3]))+...
    %         coeff(4)*(circshift(Vz,[0 -3])-circshift(Vz,[0 4]))+...
    %         coeff(5)*(circshift(Vz,[0 -4])-circshift(Vz,[0 5]))+...
    %         coeff(6)*(circshift(Vz,[0 -5])-circshift(Vz,[0 6]))+...
    %         coeff(7)*(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));
    %
    Vzx=coeff(1)*(circshift(Vz,[0 0])-circshift(Vz,[0 1]));
    for iii=2:length(coeff)
        Vzx=Vzx+coeff(iii)*(circshift(Vz,[0 -iii+1])-circshift(Vz,[0 iii]));
    end
    
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
% save('BalancedSchemeM22.mat')
save('BalancedSchemeM30.mat')