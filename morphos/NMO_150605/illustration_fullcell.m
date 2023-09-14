[~,~,~,coords,r,~]=readSWC('./refinements_ref2.swc-ref4.swc');
markerSize = (r./max(r)).*50;
% markerSize = (r./max(r)).*0 + 25;
fnt_size = 20;
LineWidth = 1.25;


figure; 
% p = scatter3(coords(:,1),coords(:,2),coords(:,3),markerSize,'filled','MarkerFaceColor','black'); 
p = scatter(coords(:,1),coords(:,2),markerSize,'filled','MarkerFaceColor','black');

p(1).LineWidth = LineWidth;
set(gca,'FontSize',fnt_size,'box','off')
[ti,s] = title('');
ti.FontSize = fnt_size;
ti = xlabel('{\mu}M');
ti.FontSize = fnt_size;
ti = ylabel('{\mu}M');
ti.FontSize = fnt_size;
ti = zlabel('{\mu}M');
ti.FontSize = fnt_size;

xlim([-40 40]); ylim([-60 80]); 
% xlim([-8.3-1.33 -8+1.33]); ylim([-10.2-1.33 -9.8+1.33]); % zlim([-4,0]);

pt1 = [-0.00022222,-0.063778,0];
pt2 = [-0.04, 20.16, -1.71];
pt3 = [0.26, 20.95, 0.28];
pt4 = [0.3, -4.63, 0.07];
pt5 = [-7.91, -9.7, -1.9];
pt6 = [-8.08, -9.91, -1.9];
pt7 = [2.17, -8.61, 0.07];
pt8 = [3.52, -11.15, 0.07];
pt9 = [3.94, -11.68, 0.07];

pt = [pt1; pt2; pt3; pt4; pt5; pt6; pt7; pt8; pt9];

hold on
% scatter3(pt1(1),pt1(2),pt1(3),70,'filled','MarkerFaceColor','r', 'MarkerEdgeAlpha',0.01)
clsch = hsv(9)
% scatter(pt1(1),pt1(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(1,:,:));

%% For 2D
scatter(pt1(1),pt1(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(1,:,:));
scatter(pt2(1),pt2(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(2,:,:));
scatter(pt3(1),pt3(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(3,:,:));
scatter(pt4(1),pt4(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(4,:,:));
scatter(pt5(1),pt5(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(5,:,:));
scatter(pt6(1),pt6(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(6,:,:));
scatter(pt7(1),pt7(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(7,:,:));
scatter(pt8(1),pt8(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(8,:,:));
scatter(pt9(1),pt9(2),70,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(9,:,:));

scatter(pt1(1),pt1(2),70,'MarkerFaceColor',clsch(1,:,:));
scatter(pt2(1),pt2(2),70,'MarkerFaceColor',clsch(2,:,:));
scatter(pt3(1),pt3(2),70,'MarkerFaceColor',clsch(3,:,:));
scatter(pt4(1),pt4(2),70,'MarkerFaceColor',clsch(4,:,:));
scatter(pt5(1),pt5(2),70,'MarkerFaceColor',clsch(5,:,:));
scatter(pt6(1),pt6(2),70,'MarkerFaceColor',clsch(6,:,:));
scatter(pt7(1),pt7(2),70,'MarkerFaceColor',clsch(7,:,:));
scatter(pt8(1),pt8(2),70,'MarkerFaceColor',clsch(8,:,:));
scatter(pt9(1),pt9(2),70,'MarkerFaceColor',clsch(9,:,:));
%% For 3D
% scatter3(pt1(1),pt1(2),pt1(3),500,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(1,:,:));
% scatter3(pt2(1),pt2(2),pt2(3),100,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(2,:,:));
% scatter3(pt3(1),pt3(2),pt3(3),100,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(3,:,:));
% scatter3(pt4(1),pt4(2),pt4(3),100,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(4,:,:));
% scatter3(pt5(1),pt5(2),pt5(3),100,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(5,:,:));
% scatter3(pt6(1),pt6(2),pt6(3),100,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(6,:,:));
% scatter3(pt7(1),pt7(2),pt7(3),100,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(7,:,:));
% scatter3(pt8(1),pt8(2),pt8(3),100,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(8,:,:));
% scatter3(pt9(1),pt9(2),pt9(3),100,'filled','MarkerFaceAlpha',0.75,'MarkerFaceColor',clsch(9,:,:));

%%
l = legend();
set( l, 'Box', 'off' )