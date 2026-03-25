load('./data/qPCR.mat');

figure(1)
clf

ABdose = repmat(flip(Locations), 3, 1);

acrBperRobPerngDNA = flip(Pooled_qPCR_Data(:, 1:3))';

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

acrBperRobPerngDNA2 = flip(Pooled_qPCR_Data(:, 4:end))';

X2 = ABdose(:);
Y2 = acrBperRobPerngDNA2(:);
grey = [1 1 1]*0.7;
x2 = 35:2:195;

errorbar(ABdose(1,:),mean(acrBperRobPerngDNA2),1.96*std(acrBperRobPerngDNA2)/sqrt(2),'linestyle','none',...
    'Color',grey,'MarkerSize',30,'MarkerFaceColor',grey)
hold on
plot(X2,Y2,'.','color',grey,'MarkerSize',12)
plot(ABdose(1,:),mean(acrBperRobPerngDNA2),'.k','MarkerSize',32)

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

