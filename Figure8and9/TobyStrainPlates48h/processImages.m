disp('This is the code that processes raw images but it will not work here ')
disp('because there is not enough repository space available to store them all. ')
disp('Example raw images can be found in the paper and in the sub-directory ./images/ ')
disp('Processed data derived from those raw images can be found in the sub-directory ./data/ ')

pause

%%

clear('variables')
close all
clc

saveFigures = 0;
plotOn = 0;

%%

rootDir = './images/';

D = dir(rootDir);
L = length({D.name});
files = {};

for l = 1:L
    n = D(l).name;
    if length(n) > 2
        files = union(files,[rootDir,n]);
    end
end

%%

dummyEntry = {'dummy',1:2,1:2,1,2,NaN};
variables = {'fileLabel','roiX','roiY','centerX','centerY','drugDose'};

drugZoneData = {...
    {'NDC',1900:2550,1700:2300,330,320,0},...
    {'DZTL',1290:1810,990:1570,245,270,300},...
    {'DZTR',2490:3150,990:1650,330,320,250},...
    {'DZBL',1420:1980,2550:3100,280,280,100},...
    {'DZBR',2710:3370,2290:2880,350,295,150},...
    };

%%

filegroups = {[1 2 3 4]};
for j = 1:191
    filegroups{j+1} = 4*j+1:4*j+4;
end

%%

acrOperonHistograms = {};
acrProteinHistograms = {};
proteinDataSize = {};
operonDataSize = {};
radialCR = {};
radialYR = {};

%%

drugZone = 4;
        
for filegrouping = 1:length(filegroups)
    
    close all
    
    %%
    
    group = filegroups{filegrouping};
    imageBF = single(imread(files{group(1)}));
    imageCFP = single(imread(files{group(2)}));
    imageYFP = single(imread(files{group(3)}));
    imageRFP = single(imread(files{group(4)}));

    %%

    if not(plotOn)
        disp(['Processing file group ',num2str(filegrouping),' of ',num2str(length(filegroups))]);        
    end

    %close all
    if not(plotOn)
        disp(['Processing drug zone ',num2str(drugZone),' of ',num2str(5)]);
    end

    %%

    thisDrugZoneData = drugZoneData{drugZone};

    fileLabel = thisDrugZoneData{1};
    roiX = thisDrugZoneData{2}; 
    roiY = thisDrugZoneData{3};
    centerX = thisDrugZoneData{4};
    centerY = thisDrugZoneData{5};
    drugDose = thisDrugZoneData{6};

    if plotOn
        figure(1)
        set(1,'pos',[680         394        1053         703])
        subplot(2,2,1);
        imagesc(uint8(imageBF));

        marker = '.k';
        hold on
        plot(min(roiX),min(roiY),marker,'markersize',20);
        plot(min(roiX),max(roiY),marker,'markersize',20);
        plot(max(roiX),min(roiY),marker,'markersize',20);
        plot(max(roiX),max(roiY),marker,'markersize',20);

        plot([min(roiX) min(roiX)],[min(roiY) max(roiY)],'-k','linewidth',1);
        plot([max(roiX) max(roiX)],[min(roiY) max(roiY)],'-k','linewidth',1);
        plot([min(roiX) max(roiX)],[min(roiY) min(roiY)],'-k','linewidth',1);
        plot([min(roiX) max(roiX)],[max(roiY) max(roiY)],'-k','linewidth',1);

        title('bright field marked with ROI')

        subplot(2,2,2);
        imagesc(uint8(imageCFP));
        title('CFP')
        subplot(2,2,3);
        imagesc(uint8(imageYFP));
        title('GFP')
        subplot(2,2,4);
        imagesc(uint8(imageRFP));
        title('RFP')
    end

    %%

    BF = imgaussfilt(imageBF(roiY,roiX,:));
    CFP = imgaussfilt(imageCFP(roiY,roiX,3));
    YFP = imgaussfilt(imageYFP(roiY,roiX,2));
    RFP = imgaussfilt(imageRFP(roiY,roiX,1));

    %%

    if plotOn
        figure(2)
        set(2,'pos',[680         394        1053         703])
        subplot(2,2,1);
        imagesc(uint8(imageBF));
        title('bright field')
        subplot(2,2,2);
        imagesc(uint8(CFP));
        title('CFP')
        subplot(2,2,3);
        imagesc(uint8(YFP));
        title('GFP')
        subplot(2,2,4);
        imagesc(uint8(RFP));
        title('RFP')
    end

    %%

    CFPperRFP = CFP ./ RFP;
    YFPperRFP = YFP ./ RFP;
    YFPperCFP = YFP ./ CFP;

    if plotOn
        figure(3)
        set(3,'pos',[543         568        1329         289])

        subplot(1,3,1)
        surf(CFPperRFP,'EdgeColor','none')
        colorbar
        axis tight
        view(2)
        title('cyan per red')

        subplot(1,3,2)
        surf(YFPperRFP,'EdgeColor','none')
        colorbar
        axis tight
        view(2)
        title('green per red')

        subplot(1,3,3)
        surf(YFPperCFP,'EdgeColor','none')
        colorbar
        axis tight
        view(2)
        title('green per cyan')
    end

    backgroundROI = 1:30;

    CPR = imgaussfilt(CFPperRFP(backgroundROI,backgroundROI));
    CPR = sort(CPR(:));
    NCPR = length(CPR);
    CPR = CPR(floor(0.2*NCPR):floor(0.8*NCPR));
    CFPperRFPbg = mean(CPR);

    YPR = imgaussfilt(YFPperRFP(backgroundROI,backgroundROI));
    YPR = sort(YPR(:));
    NYPR = length(YPR);
    YPR = YPR(floor(0.2*NYPR):floor(0.8*NYPR));
    YFPperRFPbg = mean(YPR);

    %%

    nCFPperRFP = imgaussfilt(CFPperRFP/CFPperRFPbg);
    nCFPperRFP(nCFPperRFP > 4) = NaN;
    nYFPperRFP = imgaussfilt(YFPperRFP/YFPperRFPbg);
    nYFPperRFP(nYFPperRFP > 4) = NaN;

    if plotOn
        figure(4)
        set(4,'pos',[543         568        1329         289])

        subplot(1,2,1)
        surf(nCFPperRFP,'EdgeColor','none')
        hold on
        plot3(centerX,centerY,1,'.w','markersize',30);

        axis tight
        zlim([0 4])
        view(2)
        colorbar
        title('cyan per red bg normalised')

        subplot(1,2,2)
        surf(nYFPperRFP,'EdgeColor','none')
        axis tight
        zlim([0 4])
        view(2)
        colorbar
        title('green per red bg normalised')
    end

    %%

    pumpsPerCell = nYFPperRFP(:);
    genesPerCell = nCFPperRFP(:);

    minVal = 0.5;
    bumps = 2;

    pumpsPerCell = pumpsPerCell(pumpsPerCell > minVal);
    genesPerCell = genesPerCell(genesPerCell > minVal);

    [xg,hg] = hist(genesPerCell,60);
    [xp,hp] = hist(pumpsPerCell,60);

    if plotOn
        figure(5)
        set(gcf,'pos',[543         568        1329         289])

        subplot(1,2,1)
        bar(hg,xg/sum(xg))
        xlabel('acr operons per cell proxy')
        ylabel('frequency')
        xlim([minVal 3])
        set(gca,'Xtick',[minVal 1 2 3])

        options = statset('Display','final');
        gm = fitgmdist(genesPerCell,bumps,'Options',options);
        gmPDF = @(x)pdf(gm,x);
        pltPDF = gmPDF(hg')';
        pltPDF = pltPDF/sum(pltPDF);
        hold on
        plot(hg,pltPDF);
        legend({['image data (n = ',num2str(length(genesPerCell)),')'],'Gaussian mixture model'});

        subplot(1,2,2)
        bar(hp,xp/sum(xp));
        xlabel('AcrB per cell proxy')
        ylabel('frequency')
        xlim([minVal 3])
        set(gca,'Xtick',[minVal 1 2 3])

        options = statset('Display','final');
        gm = fitgmdist(pumpsPerCell,bumps,'Options',options);
        gmPDF = @(x)pdf(gm,x);
        pltPDF = gmPDF(hp')';
        pltPDF = pltPDF/sum(pltPDF);
        hold on
        plot(hp,pltPDF);
        legend({['image data (n = ',num2str(length(pumpsPerCell)),')'],'2-Gaussian model'});
        drawnow
    end

    acrProteinHistograms{filegrouping} = [hp ; xp/sum(xp)];
    proteinDataSize{filegrouping} = length(genesPerCell);

    acrOperonHistograms{filegrouping} = [hg ; xg/sum(xg)];
    operonDataSize{filegrouping} = length(pumpsPerCell);

    zoneLabel{filegrouping} = fileLabel;
    dose{filegrouping} = drugDose;

    %%

    N = 500;
    H = 100;

    theta = 2*pi*((-H):(H-1))./H;
    R = 130:N;

    dataCR = zeros(size(R));
    dataYR = zeros(size(R));

    if plotOn
        figure(4)
        subplot(1,2,1)
        xL = xlim;
        yL = ylim;
        zL = zlim;

        hold on
    end

    for k = 1:length(R)
        r = R(k);

        datumCR = 0*theta;
        datumYR = 0*theta;

        for j = 1:length(theta)
            th = theta(j);
            posX = uint16(centerX + r*cos(th));
            posY = uint16(centerY + r*sin(th));
            try
                datumCR(j) = nCFPperRFP(posX,posY);
            catch
                datumCR(j) = NaN;
            end
            try
                datumYR(j) = nYFPperRFP(posX,posY);
            catch
                datumYR(j) = NaN;
            end
        end

        datumCR = datumCR(~isnan(datumCR));
        dataCR(k) = mean(datumCR);

        datumYR = datumYR(~isnan(datumYR));
        dataYR(k) = mean(datumYR);

        try
            X = uint16(centerX + r*cos(theta));
            Y = uint16(centerY + r*sin(theta));
            if plotOn
                plot3(X,Y,2*ones(size(X)),'.w','markersize',4)
            end
        catch
            disp('error line 335')
        end

    end

    for k = 1:3
        s = k / 10;
        thisr = 1 + R(round(s*N)) - R(1);
        X = uint16(centerX + thisr*0);
        Y = uint16(centerY + thisr*1);
        try
            if plotOn
                plot3(X,Y,3,'.k','markersize',20)
            end
            text(double(1.025*X),double(1.025*Y),3.5,num2str(thisr));
        end
    end

    if plotOn
        xlim(xL)
        ylim(yL)
        zlim(zL)
    end

    %%

    if plotOn
        figure(6)
        plot(R,dataCR)
        hold on
        plot(R,dataYR)
        axis tight
        legend({'acr operons per cell proxy','AcrB per cell proxy'});
        xlabel('distance from drug centre')
        ylabel('units')
        drawnow
    end

    radialCR{filegrouping} = dataCR;
    radialYR{filegrouping} = dataYR;

    %%

    if plotOn
        if saveFigures
            for f = 1:6
                figure(f)
                export_fig(['analysisFigure-',fileLabel,'-',num2str(f),'.jpg']);
            end
        end
    end

end

%%

save([fileLabel,'dataset.mat'],...
    'zoneLabel','dose',...
    'acrProteinHistograms','acrOperonHistograms',...
    'operonDataSize','proteinDataSize','radialCR','radialYR','R')

