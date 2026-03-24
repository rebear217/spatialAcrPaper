clear variables
clc
close all

%%

files = {'NDCdataset','DZTRdataset','DZTLdataset','DZBLdataset','DZBRdataset'};

Xmins = [0.65,0.65,0.65,0.65,0.9];
Xmaxes = [1.9,2.5,2.5,2.5,2.5];
saveFigures = 1;

times = 1:191;
Dose = zeros(1,5);
for f = 1:5
        
    load(['./data/',files{f},'.mat']);    
    Dose(f) = dose{1};
    
    Xmax = Xmaxes(f);
    Xmin = Xmins(f);    
    x{f} = Xmin:0.0125:Xmax;    
    
    R = radialCR{1};
    R = R(~isnan(R));
    
    space{f} = 1:length(R);
    
    AcrBHistoSurface{f} = zeros(length(times),length(x{f}));
    acrHistoSurface{f} = zeros(length(times),length(x{f}));

    meansOperon{f} = zeros(length(times),1);
    stdsOperon{f} = zeros(length(times),1);

    meansProtein{f} = zeros(length(times),1);
    stdsProtein{f} = zeros(length(times),1);

    for n = times
        X = acrProteinHistograms{n}(1,:);
        Y = acrProteinHistograms{n}(2,:);

        H = interp1(X,Y,x{f},'spline','extrap');    
        AcrBHistoSurface{f}(n,:) = H.*(H>0);
        
        m = sum(X.*Y);
        v = sum((X-m).^2.*Y);
        s = sqrt(v);
        meansProtein{f}(n) = m;
        stdsProtein{f}(n) = v;
        
        X = acrOperonHistograms{n}(1,:);
        Y = acrOperonHistograms{n}(2,:);

        H = interp1(X,Y,x{f},'spline','extrap');    
        acrHistoSurface{f}(n,:) = H.*(H>0);
        
        m = sum(X.*Y);
        v = sum((X-m).^2.*Y);
        s = sqrt(v);
        meansOperon{f}(n) = m;
        stdsOperon{f}(n) = s;
    end

end

%%

[dose,DI] = sort(Dose);

for j = 1:2
    f = DI(j);
    figure()
    set(gcf,'pos',[986         647        700         567])
    if j == 1
        S = 0:0.1:3;
        T = 0:0.1:47.7;
        surf(S,T,zeros(length(T),length(S)),'EdgeColor','none')
        hold on
    end
    surf(x{f},times/4,AcrBHistoSurface{f},'EdgeColor','none')
    view(2)
    axis tight
    xlim([0.7 2.5])
    zlim([0.0 0.3])
    
    ylabel('time (h)','FontSize',18)
    xlabel('AcrB per genome','FontSize',18)
    legend(['dose ',num2str(Dose(f)),'x MIC'],'Color','w','FontSize',18);
    colorbar('location','EastOutside');
    colormap(hot.^(1))
end


%%

ss = 10;
lw = 4;

for j = 1:2
    f = DI(j);
    figure()
    set(gcf,'pos',[986         647        700         567])
    L = length(AcrBHistoSurface{f}(1,:));
    for n = 10:floor(L/15):L
        s = (n-1)/(L-1);
        c = [s/2 0.2 s];
        X = n*ones(size(x{f}))/4;
        Y = x{f};
        Z = smooth(AcrBHistoSurface{f}(n,:),ss);
        plot3(X,Y,Z,'color',c,'linewidth',lw);
        hold on
    end
	axis tight
    %ylim([0.5 2.2])
    %zlim([0.0 0.3])
    view(-210,32.5)
    
    xlabel('time (h)','FontSize',18)
    ylabel('AcrB per genome','FontSize',18)
    zlabel('frequency','FontSize',18)
    
    title(['dose ',num2str(Dose(f)),'xMIC']);
    box on
end

%%

figure()
set(gcf,'pos',[296         664        1681         452])

lw = 2;
ss = 10;

legs = {'mean acr per genome @ 0xMIC','100xMIC','150xMIC','250xMIC','300xMIC'};
legs2 = {'mean AcrB per genome @ 0xMIC','100xMIC','150xMIC','250xMIC','300xMIC'};
legs3 = {'mean AcrB per acr @ 0xMIC','100xMIC','150xMIC','250xMIC','300xMIC'};

for f = 1:5
    cf = num2str(corr(meansProtein{f},meansOperon{f}),3);
    legs3{f} = [legs3{f},' (corr coeff ',cf,')'];
end

for j = 1:5

    f = DI(j);
    s = (j-1)/4;
    c = [s s/4 1-s];
    
    subplot(1,3,1)
    plot(times/4,smooth(meansOperon{f},ss),'color',c,'linewidth',lw,'DisplayName',legs{j});
    hold on

    T = times/4;
    Y = smooth(meansOperon{f},ss);
    E = 100*stdsOperon{f}/sqrt(operonDataSize{f}-1);
    J = 10;
    errorbar(T(1:J:end),Y(1:J:end),E(1:J:end),'.','color',c,'linewidth',1,'HandleVisibility','off');

    xlabel('time (h)');
    ylabel('image histogram means')
    axis tight
    ylim([0.8 2])
    xL = xlim;
    xlim([1 25])
    legend boxoff
    legend('location','northwest')

    subplot(1,3,2)
    T = times/4;
    Y = smooth(meansProtein{f},ss);
    E = 100*stdsProtein{f}/sqrt(proteinDataSize{f}-1);
    J = 10;
    plot(times/4,smooth(meansProtein{f},ss),'-','color',c,'linewidth',lw,'DisplayName',legs2{j});
    hold on
    errorbar(T(1:J:end),Y(1:J:end),E(1:J:end),'.','color',c,'linewidth',1,'HandleVisibility','off');
    xlabel('time (h)');
    ylabel('image histogram means')
    axis tight
    ylim([0.8 2])
    xL = xlim;
    xlim([1 25])
    %legend('mean AcrB per genome @ 0xMIC','100xMIC','150xMIC','250xMIC','300xMIC')
    legend boxoff
    legend('location','northwest')
    
    subplot(1,3,3)
    plot(times/4,smooth(meansProtein{f}./meansOperon{f},ss),'color',c,'linewidth',2,'DisplayName',legs3{j});
    hold on
    xlabel('time (h)');
    ylabel('image histogram means')
    axis tight
    ylim([0.8 2])
    xL = xlim;
    xlim([1 25])
    legend boxoff
    legend('location','northwest')
end

%%

figure()
set(gcf,'pos',[765   665   467   381])

N1 = 100;
N2 = 191;

for j = N1:1:N2

    s = (j-N1)/(N2-N1);
    c = [s s/4 1-s];
    
    myVector = [0,0,0,0,0];
    for k = 1:5
        f = DI(k);
        v = meansOperon{f};
        myVector(k) = v(j);
    end
    
    if j > N1 && j < N2
        plot(dose,myVector,'-','color',c,'linewidth',1,'HandleVisibility','off');
    else
        plot(dose,myVector,'-','color',c,'linewidth',1,'DisplayName',['times closer to ',num2str(T(j)),'h']);
    end

    hold on
    drawnow
end

xlabel('dose (xMIC)');
ylabel('acr image histogram means')
set(gca,'Xtick',dose)

axis tight
ylim([0.9 1.8])
xL = xlim;
xlim([0 max(dose)])
legend()
legend boxoff
legend('location','northwest')    

%exportgraphics(gcf,'./figures/meanDoseResponse.pdf')

%%

figure()
set(gcf,'pos',[986   657   604   557])

c = {'k','b'};
minTs = [14,14];

for j = 1:2
    f = DI(j);

    Xs = [];
    Ys = [];

    L = length(AcrBHistoSurface{f}(1,:));
    for n = 1:1:L
        X = n*ones(size(x{f}))/4;
        Y = x{f};
        Z = smooth(AcrBHistoSurface{f}(n,:),10);
        zm = mean(Z);
        if mod(n,10) == -1
            plot3(X,Y,Z,'color',c{j},'linewidth',2);
        end
        iLM = islocalmax(Z) & (Z>1.75*zm);
        zLM = Z(iLM);
        xLM = X(iLM);
        yLM = Y(iLM);
        [zLM,I] = sort(zLM);
        xLM = xLM(I);
        yLM = yLM(I);
        N = min(2,length(xLM));

        if min(xLM) > minTs(j)
            ms = 30;
            sp = '.';
        else
            ms = 2;
            sp = 'o';
        end
        if n == L
            plot(xLM(1:N),yLM(1:N),sp,'markersize',ms,'color',c{j},'DisplayName',['local maxima @ dose ',num2str(Dose(f)),'xMIC'])
        else
            plot(xLM(1:N),yLM(1:N),sp,'markersize',ms,'color',c{j},'HandleVisibility','off')
            hold on
        end

        Xs = [Xs xLM(1:N)];
        Ys = [Ys yLM(1:N)];        
    end
	axis tight
    ylim([0.7 2])
    %zlim([0.0 0.3])
    %view(-150,33)
    
    xlabel('time (h)')
    ylabel('AcrB per genome local maxima')
    %zlabel('frequency')
    mYs = movmean(Ys,15);
    plot(Xs(Xs < minTs(2)),mYs(Xs < minTs(2)),'.','markersize',30,'color',c{j},'HandleVisibility','off');

end

legend('location','northwest')
xlim([0 30])

%%


figure()
set(gcf,'pos',[986   657   604   450])

j = 2;
f = DI(j);

Xs = [];
Ys = [];
Zs = [];

Ns = 11;

L = length(AcrBHistoSurface{f}(1,:));
for n = [Ns:25 110:L]

    s = (n-1)/L;
    c = [s/2 0.2 s];      
    
    lw = 1 + 1*s;
    if n <= Ns || n == L
        lw = 6;
    end

    X = n*ones(size(x{f}))/4;
    Y = x{f};
    Z = smooth(AcrBHistoSurface{f}(n,:),10);
    zm = mean(Z);

    if n == Ns
        plot(Y,Z,'color',c,'linewidth',lw,'DisplayName','early distributions');
    elseif n == L
        plot(Y,Z,'color',c,'linewidth',lw,'DisplayName','late distributions');
    else
        plot(Y,Z,'color',c,'linewidth',lw,'HandleVisibility','off');
    end
    hold on       

end
axis tight
%ylim([0.7 2])
%xlim([0 30])

ylabel('frequency')
xlabel('AcrB per genome')

legend('location','northeast')


%%

if saveFigures
    for fg = 1:8
        figure(fg)
        export_fig(['./figures/Figure',num2str(fg),'.jpg']);
    end
end