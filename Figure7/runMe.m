figure(1)
clf

ABdose = [40 130 160 190;...
     40 130 160 190;...
     40 130 160 190];

acrBperRobPerngDNA = [2.63 1.49 1.71 1.06;...
                      1.79 1.32 1.02 0.97;...
                      1.74 1.28 1.05 0.91];

X = ABdose(:);
Y = acrBperRobPerngDNA(:);
grey = [1 1 1]*0.7;
x = [35 195];
errorbar(ABdose(1,:),mean(acrBperRobPerngDNA),1.96*std(acrBperRobPerngDNA)/sqrt(2),'linestyle','none',...
    'Color',grey,'MarkerSize',30,'MarkerFaceColor',grey)
hold on
plot(X,Y,'.','color',grey,'MarkerSize',12)
plot(ABdose(1,:),mean(acrBperRobPerngDNA),'.k','MarkerSize',32)

axis tight
plot([35 190],[1 1],'-k')
ylim([0 3.1])
xlim(x)
xlabel('Distance to Antbiotic (px)')
ylabel('${\sf acrB} \cdot {\sf rob}^{-1} \cdot {\sf ng}^{-1} {\sf DNA}$','Interpreter','latex')

F = fitlm(X,Y);
plot(x,F.feval(x),'-r')

%%

figure(2)
clf

ABdose2 = [40 70 130 160 190;...
     40 70 130 160 190;...
     40 70 130 160 190];

acrBperRobPerngDNA2 = [1.62 3.04 1.81 1.08 1.07;...
                      1.47 1.81 1.88 0.91 1.03;...
                      1.35 2.4 1.82 0.61 0.81];

X2 = ABdose2(:);
Y2 = acrBperRobPerngDNA2(:);
grey = [1 1 1]*0.7;
x2 = 35:2:195;

errorbar(ABdose2(1,:),mean(acrBperRobPerngDNA2),1.96*std(acrBperRobPerngDNA2)/sqrt(2),'linestyle','none',...
    'Color',grey,'MarkerSize',30,'MarkerFaceColor',grey)
hold on
plot(X2,Y2,'.','color',grey,'MarkerSize',12)
plot(ABdose2(1,:),mean(acrBperRobPerngDNA2),'.k','MarkerSize',32)

axis tight
plot([35 190],[1 1],'-k')
ylim([0 3.1])
xlim(x)
xlabel('Distance to Antbiotic (px)')
ylabel('${\sf acrB} \cdot {\sf rob}^{-1} \cdot {\sf ng}^{-1} {\sf DNA}$','Interpreter','latex')

F2 = fitlm(X2,Y2);
plot(x,F2.feval(x),'-','color',grey)

Q = fitnlm(X2,Y2,@(b,x)(b(1)+b(2)*x+b(3)*x.^2),[1 1 1])
plot(x2,Q.feval(x2),'-r')

