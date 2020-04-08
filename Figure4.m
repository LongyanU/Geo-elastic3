
load('BalancedSmallGridBoundary.mat')
Font_size = 20;
xlabel('x/dx','FontSize',Font_size)
ylabel('z/dz','FontSize',Font_size)

imagesc(vp(46*2:end,:))
hold on; plot(200*2,100*2,'r*','linewidth',2)

hold on; plot(250*2,55*2,'bo','linewidth',2)
hold on; plot(300*2,55*2,'bo','linewidth',2)
hold on; plot(400*2,55*2,'bo','linewidth',2)
hold on; plot(500*2,55*2,'bo','linewidth',2)
hold on; plot(600*2,55*2,'bo','linewidth',2)
hold on; plot(650*2,55*2,'bo','linewidth',2)

hold on; plot(400*2,(150-46)*2,'bo','linewidth',2)
hold on; plot(400*2,(200-46)*2,'bo','linewidth',2)
hold on; plot(400*2,(250-46)*2,'bo','linewidth',2)