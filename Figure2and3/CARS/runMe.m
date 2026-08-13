clear all
close all
clc

%% plots Figure 2C (all subplots) in the spatialAcr paper:

close all
for T = [10,15,20]
    outputT = timesolvePanel(T);
    close all
end

%% plots Figure 2A-B in the spatialAcr paper:

for option = 1
    close all
    [mySurface,parameters,mySSurface,myRSurface] = runTheModel(option);
    close all
    plotWave(mySurface,parameters,num2str(option));
    plotFinalWave(mySSurface,myRSurface,parameters,num2str(option))
end

%% plots Figure 3A-B in the spatialAcr paper:

for option = 2:8
    close all
    [mySurface,parameters,mySSurface,myRSurface] = runTheModel(option);
    close all
    plotFinalWave(mySSurface,myRSurface,parameters,num2str(option))
end

