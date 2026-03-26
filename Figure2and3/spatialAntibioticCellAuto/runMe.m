close all
clear all
clc

%%

PLOTSAVEON = 1;
Tmax = 70;
[state,plotStates] = simulation(Tmax,PLOTSAVEON);

%%

figure(1)
clf
s=4;

surf(plotStates(s:end,:),'EdgeColor','none','FaceColor','interp')
axis tight
box on
camlight right
lighting gouraud
alpha(1)

set(gca,'Xtick',0:20:100)
set(gca,'XtickLabel',100:-20:0)

ylabel('time (model epochs)')
xlabel('distance from drug source')
zlabel('bacterial population density')

%lighting gouraud
X = 1:100;
hold on
for T = 1:10:Tmax
    if T == 1
        plot3(X,T*ones(size(X)),plotStates(T+s-1,:),'-','Color','r',...
            'LineWidth',3)
    else
        plot3(X,T*ones(size(X)),plotStates(T+s-1,:),'-','Color',[1 1 1]*0.8,...
            'LineWidth',3)
    end
end
%colorbar
view(6.25,53.4)

figure(2)
clf
for T = 1:10:Tmax
    S = (T-1)/Tmax;
    c = 0.7*[1 1 1]*(1-S) + [0 0 0]*S;
    lw = 3*S + 1;
    if T == 1
        plot(fliplr(X),plotStates(T+s-1,:),'-','Color',[1 0 0],'LineWidth',2,...
            'DisplayName',['pop density profile @T=',num2str(T),' model epoch'])
    else
        plot(fliplr(X),plotStates(T+s-1,:),'-','Color',c,'LineWidth',lw,...
            'DisplayName',[num2str(T),' epochs'])
    end
    hold on
end
legend('Location','northwest')

xlabel('distance from drug source')
ylabel('bacterial population density')

%%

figure(1)
exportgraphics(gcf,'./figures/CAwave3d.pdf')
figure(2)
exportgraphics(gcf,'./figures/CAwave2d.pdf')

