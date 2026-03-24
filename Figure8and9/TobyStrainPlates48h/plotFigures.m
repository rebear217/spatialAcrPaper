clear variables
clc
close all

%%

files = {'NDCdataset','DZTRdataset','DZTLdataset','DZBLdataset','DZBRdataset'};
sFactor = 50;

Tmin = 5;

%%

for f = 1:5
    
    close all
    
    load(['./data/',files{f},'.mat']);

    Xmins = [0.65,0.65,0.65,0.65,0.9];
    Xmaxes = [1.5,2.5,2.5,2.5,2.5];
    Xmax = Xmaxes(f);
    Xmin = Xmins(f);

    saveFigures = 1;

    times = 1:100;

    R = radialCR{1};
    R = R(~isnan(R));
    space = 1:205;

    CRsurface = zeros(length(times),length(space));
    YRsurface = zeros(length(times),length(space));

    for n = times
        %CRsurface(n,:) = flipud(smooth(radialCR{n}(space),'sgolay'));
        %YRsurface(n,:) = flipud(smooth(radialYR{n}(space),'sgolay'));

        CRsurface(n,:) = radialCR{n}(space);
        YRsurface(n,:) = radialYR{n}(space);
    end

    for r = space
        CRsurface(:,r) = smooth(CRsurface(:,r));
        YRsurface(:,r) = smooth(YRsurface(:,r));
    end

    %%

    figure(1)
    set(1,'pos',[786         270        1477         960])

    subplot(2,2,1)

    surf(times/4,space,CRsurface','EdgeColor','none')
    view(2)
    axis tight
    colorbar
    zlim([0.8 2.5])
    ylabel('distance from antibiotic (px)')
    xlabel('time (h)')
    title('acr operons per genome proxy')
    colormap(hot)
    xL = xlim;
    xlim([4 xL(2)])

    subplot(2,2,2)
    surf(times/4,space,YRsurface','EdgeColor','none')
    view(2)
    axis tight
    zlim([0.8 2.5])
    colorbar
    ylabel('distance from antibiotic (px)')
    xlabel('time (h)')
    title(['AcrB per genome proxy (dose ',num2str(dose{1}),'x MIC)'])
    xL = xlim;
    xlim([4 xL(2)])

    subplot(2,2,3)
    spaceslice = space(50:50:200);
    for j = 1:length(spaceslice)
        s = spaceslice(j);
        S = (j-1) / (length(spaceslice)-1);
        %c = [1-S/2,0,S];
        c = [S S/3 1-S];
        cr = mean(CRsurface(:,s-5:s+5),2);
        %cr = smooth(cr);
        plot(times/4,cr,'color',c)
        hold on
    end
    xlabel('time (h)')
    ylabel('acr operons per genome proxy')
    axis tight
    xL = xlim;
    xlim([Tmin xL(2)])
    ylim([0.8 2.5])
    
    L2 = {};    
    L = num2str(spaceslice(:));
    for j = 1:length(spaceslice)
        L2{j} = [L(j,:),' px'];
    end
    L2{1} = [L2{1},' from antibiotic'];
    legend(L2)
    legend('location','northwest')
    
    subplot(2,2,4)
    spaceslice = space(50:50:200);
    for j = 1:length(spaceslice)
        s = spaceslice(j);
        S = (j-1) / (length(spaceslice)-1);
        %c = [1-S/2,0,S];
        c = [S S/3 1-S];
        yr = mean(YRsurface(:,s-5:s+5),2);
        %yr = smooth(yr);
        plot(times/4,yr,'color',c)
        hold on
    end
    xlabel('time (h)')
    ylabel('AcrB per genome proxy')
    axis tight
    xL = xlim;
    xlim([Tmin xL(2)])
    ylim([0.8 2.5])
    
    legend(L2)    
    legend('location','northwest')
    
    %%

    x = Xmin:0.0125:Xmax;

    AcrBHistoSurface = zeros(length(times),length(x));
    acrHistoSurface = zeros(length(times),length(x));

    meansOperon = zeros(length(times),1);
    stdsOperon = zeros(length(times),1);

    meansProtein = zeros(length(times),1);
    stdsProtein = zeros(length(times),1);

    for n = times
        X = acrProteinHistograms{n}(1,:);
        Y = acrProteinHistograms{n}(2,:);

        H = interp1(X,Y,x,'pchip','extrap');    
        AcrBHistoSurface(n,:) = H.*(H > 0);    

        m = sum(X.*Y);
        v = sum((X-m).^2.*Y);
        s = sqrt(v);% / sqrt(proteinDataSize{n}-1);
        meansProtein(n) = m;
        stdsProtein(n) = s;    

        X = acrOperonHistograms{n}(1,:);
        Y = acrOperonHistograms{n}(2,:);

        H = interp1(X,Y,x,'pchip','extrap');    
        acrHistoSurface(n,:) = H.*(H > 0);    

        m = sum(X.*Y);
        v = sum((X-m).^2.*Y);
        s = sqrt(v);% / sqrt(operonDataSize{n}-1);
        meansOperon(n) = m;
        stdsOperon(n) = s;
    end

    %% 

    figure(2)
    set(2,'pos',[1000         892        1331         446])

    subplot(1,2,1)
    surf(x,times/4,AcrBHistoSurface,'EdgeColor','none')
    colorbar
    xlabel('AcrB per genome proxy')
    title('frequency')
    ylabel('time (h)')
    axis tight
    view(2)
    colormap(hot.^(1/2))

    subplot(1,2,2)
    surf(x,times/4,acrHistoSurface,'EdgeColor','none')
    colorbar
    xlabel('acr operons per genome proxy')
    title(['frequency (dose ',num2str(dose{1}),'x MIC)'])
    ylabel('time (h)')
    axis tight
    view(2)
    colormap(hot.^(1/2))

    %%

    figure(3)
    set(3,'pos',[1021         436        1371         527])

    smP = meansProtein;
    smO = meansOperon;

    subplot(1,2,1)
    p=plot(times/4,smP);
    hold on
    q=plot(times/4,smO);
    errorbar(times(1:10:end)/4,smP(1:10:end),stdsProtein(1:10:end),'ok','Markersize',1,'color',p.Color);
    errorbar(times(1:10:end)/4,smO(1:10:end),stdsOperon(1:10:end),'ok','Markersize',1,'color',q.Color);
    
    xlabel('time (h)')
    ylabel('image histogram means')
    legend('AcrB proxy \pm std','acr proxy \pm std','location','northwest')
    axis tight
    ylim([0.8 3])

    errPO = abs(stdsProtein ./ smO) + abs(smP .* stdsOperon ./ (smO.^2));

    subplot(1,2,2)
    plot(times/4,smP ./ smO)
    hold on
    errorbar(times(1:10:end)/4,smP(1:10:end) ./ smO(1:10:end),errPO(1:10:end),'ok','Markersize',1,'color',p.Color);

    xlabel('time (h)')
    ylabel('image histogram means')
    legend(['AcrB per acr proxy \pm std (dose ',num2str(dose{1}),'x MIC)'],'location','northwest')
    axis tight
    ylim([0.8 3])

    %%

    if saveFigures
        for fg = 1:3
            figure(fg)
            export_fig(['./figures/',files{f},'Figure-',num2str(fg),'.jpg']);
        end
    end

end

%%

