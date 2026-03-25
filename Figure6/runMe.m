%% Produce Figure 6A
figure(1);
for i = 1:2
  if i == 2
    load('./data/6A_8h.mat');
  else
    load('./data/6A_50h.mat');
  end

  % Compute and plot wavefronts
  MBC_SE = std(MBC) ./ sqrt(Replicates - 1);
  MBCs(i, :) = MBC;
  
  line([mean(MBC) mean(MBC)], [0 max(max(OD_Data)) * 0.1], 'Color', 'k');
  hold on;
  line([mean(MBC) - MBC_SE*1.96, mean(MBC) + MBC_SE*1.96],...
       [max(max(OD_Data)) * 0.1, max(max(OD_Data)) * 0.1],...
       'LineWidth', 4, 'Color', 'k');
  % Highligh ring 1 at t=8h
  if i == 1
      smooth_Color = [0.7, 0.9 1] * 0.9;
      p = patch([155; 155; 175; 175], [0; 1.2; 1.2; 0], [1 1 1] * 0.95);
          set(p, 'LineStyle', 'none');
      annotation('doublearrow', [0.745 0.745], [0.65 0.1495], 'Head1Style', 'ellipse',...
                 'Head1Length', 1, 'Head1Width', 1, 'Head2Length', 4, 'Head2Width', 4,...
                 'Color', 'k');
      annotation('arrow', [0.82 0.745], [0.81 0.81], 'HeadLength', 4,...
                 'HeadWidth', 4, 'Color', 'k');
      text(165, 0.35, 'Growth Ring 1', 'Rotation', 90, 'Color', 'k',...
           'FontSize', 14, 'HorizontalAlignment', 'center');
      text(165, 1.065, ['More', sprintf('\n'), 'drug'],...
          'HorizontalAlignment', 'center', 'FontSize', 12);
  else
      smooth_Color = [0.3 0.75 0.9] * 0.9;
  end

  % Plot data
  sm = plot(Px_OD, smooth(mean(OD_Data), 11, 'sgolay'), 'LineWidth', 2,...
            'Color', smooth_Color);
  od = plot(Px_OD, mean(OD_Data), 'LineWidth', 1, 'Color', [1 1 1] * 0);

  % 95% CI
  plot(Px_OD, mean(OD_Data) - OD_SE * 1.96, Px_OD, mean(OD_Data) + OD_SE * 1.96,...
       'LineWidth', 0.5, 'Color', [1 1 1] * 0.35, 'LineStyle', ':');
    
  % Annotations
  if i == 1
      text(118, 0.1, ['WFT @ 8h:', sprintf('\n'),...
                    sprintf('%0.2f', mean(MBC)), ' \pm ',...
                    sprintf('%0.2f', MBC_SE*1.96), ' px'],...
                    'HorizontalAlignment', 'right')
  else
      text(65, 0.1, ['WFT @ 50h:', sprintf('\n'),...
                    sprintf('%0.2f', mean(MBC)), ' \pm ',...
                    sprintf('%0.2f', MBC_SE*1.96), ' px'],...
                    'HorizontalAlignment', 'right');
  end
end

% Axes
axis tight; box on;
ylim([0, 1.15]);
% Annotations
annotation(fig_3, 'doublearrow',[0.4315 0.65], [0.4 0.4],... 
           'Head2Length', 1, 'Head2Width', 1, 'Head2Style', 'ellipse',...
           'Head1Length', 4, 'Head1Width', 4);
text(97.5, 0.445, ['\Delta 42h']);
% Legend
l = legend([od; sm], {['Mean density \pm 95% CI']; 'Filtered Data'});
    set(l, 'Location', 'NorthWest', 'box', 'off', 'FontSize', 12);
ylabel(['Normalised culture density']);
xlabel('Distance from antibiotic source (px)');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Layer', 'top');
hold off;
export_fig(2, './figures/6b.pdf');


%% Produce Figure 6B
figure(2);
load('./data/6B.mat');

% Highlight standard incubation time
p0 = patch([16; 16; 24; 24], [min(min(MBCs(3:40, :))) * 0.95; max(max(MBCs(3:40, :)));...
              max(max(MBCs(3:40, :))); min(min(MBCs(3:40, :))) * 0.95], [1 1 1] * 0.9);
        set(p0, 'LineStyle', 'none');
hold on; box on;

% Plot data
p1 = plot(t(1, 3:40), mean(MBCs(3:40, :)'), 'LineWidth', 2, 'Color', 'k');
% 95% CI
plot(t(1, 3:40), mean(MBCs(3:40, :)') - MBCs_SE * 1.96,...
     t(1, 3:40), mean(MBCs(3:40, :)') + MBCs_SE * 1.96, 'Color', [1 1 1] * 0.65);
axis tight;

% Highlight wavefront at t=16h and t=24h
line([t(1, 9), t(1, end)], [mean(MBCs(9, :)), mean(MBCs(9, :))],...
      'LineStyle', ':', 'Color', 'k');

line([t(1, 13), t(1, end)], [mean(MBCs(13, :)), mean(MBCs(13, :))],...
      'LineStyle', ':', 'Color', 'k');

% Axes
xlabel('Time (h)');
ylabel('Distance to antibiotic source (px)');
set(gca, 'YDir', 'Reverse');

% Annotatios
text(5.5, 186.5, ['Edge plate'], 'FontSize', 10);
annotation('arrow', [0.21, 0.14], [0.14, 0.14], 'HeadStyle', 'rectangle',...
           'HeadLength', 0.2, 'HeadWidth', 8);
% Legend
l = legend([p1;p0], {'Mean distance \pm 95% CI'; 'Standard Incubation Time'});
        set(l, 'Location', 'SouthEast', 'box', 'off', 'FontSize', 12);
    hold off;
export_fig(2, './figures/6b.pdf');

%% Produce Figure 6C
load('./data/6C.mat');
figure(3);
l_points = [];
for i = 13:6:37
    if i == 13
        col = [0.8 0 0.2];
    elseif i == 37
        col = [1 1 1] * 0.55;
    else
        col = [1 1 1] * 0.9;
    end
    % Plot 95% CI
    plot(dist(35:end, 1), rRFUp_Mean(35:end, i) + rRFUp_SE(35:end, i) * 1.96,...
     'LineWidth', 0.75, 'LineStyle', ':', 'Color', col * 0.45);
    hold on;
    plot(dist(35:end, 1), rRFUp_Mean(35:end, i) - rRFUp_SE(35:end, i) * 1.96,...
     'LineWidth', 0.75, 'LineStyle', ':', 'Color', col * 0.45);
    % Plot smooth data
    p1 = plot(dist(35:end, 1), smooth(rRFUp_Mean(35:end, i), 9, 'sgolay'), 'LineWidth', 3,...
         'Color', col);
    % Plot raw data
    plot(dist(35:end, 1), rRFUp_Mean(35:end, i), 'LineWidth', 1, 'Color', col * 0.75);
    if i == 13 | i == 37
        set(p1, 'DisplayName', ['', num2str(t(1, i)), 'h, mean \pm 95% CI'])
        l_points = [l_points; p1];
    end
end

% Highlight data at t=24h
i = 13;
plot(dist(35:end, 1), rRFUp_Mean(35:end, i) + rRFUp_SE(35:end, i) * 1.96,...
     'LineWidth', 0.75, 'LineStyle', ':', 'Color', [0.8 0 0.2] * 0.45);
    hold on;
    plot(dist(35:end, 1), rRFUp_Mean(35:end, i) - rRFUp_SE(35:end, i) * 1.96,...
     'LineWidth', 0.75, 'LineStyle', ':', 'Color', [0.8 0 0.2] * 0.45);
    p1 = plot(dist(35:end, 1), smooth(rRFUp_Mean(35:end, i), 9, 'sgolay'), 'LineWidth', 3,...
         'Color', [0.8 0 0.2]);
    plot(dist(35:end, 1), rRFUp_Mean(35:end, i), 'LineWidth', 1, 'Color', [0.8 0 0.2] * 0.75);

% Axes
box on; set(gca, 'Layer', 'top');
axis tight; hold off;
xlabel('Distance to Antibiotic (px)');
ylabel('Relative{\it gfp-AcrB} per cell');
ylim([0, 5]);

% Legend
l = legend(l_points);
l_pos = get(l, 'Position');
    set(l, 'Location', 'NorthEast', 'box', 'off');
    set(l, 'Position', [l_pos(1)* 0.9, l_pos(2:end)]);
% Annotations
annotation('arrow', [0.4 0.6], [0.165 0.165], 'HeadLength', 6,...
           'HeadWidth', 6, 'HeadStyle', 'plain', 'LineWidth', 1);
text(110, 0.85, ['{\it acr} duplication', sprintf('\n'), 'wave back'],...
     'FontSize', 16, 'HorizontalAlignment', 'center');
export_fig(3, './figures/6c.pdf');


%% Produce Figure 6D
figure(4);
load('./data/6D.mat');

% Plot data
surf(t(1, 1:40), dist(:, 1)', rRFUp_Mean, 'LineStyle', 'none');
view(2); hold on; axis tight;
colormap(hot);

% Axes
xlabel('Time (h)');
ylabel('Distance to antibiotic (px)');
set(gca, 'TickDir', 'out');

% Colorbar
colorbar([0.925, 0.188, 0.02, 0.674]);
caxis([0.75, 3.25]);
export_fig(4, './figures/6d.pdf');
