function [R] = plotGeometry(d3bladeX,d3bladeY,d3bladeZ,d3wakeX,d3wakeY,d3wakeZ,d3wakeG,dhubHeight)

figure('Name','Geometry','NumberTitle','off')
clf
% set(gcf, 'Position', get(0, 'Screensize'));
hold on; grid on;
for i = 1:size(d3bladeX,3)
    surf(d3bladeX(:,:,i),d3bladeY(:,:,i),d3bladeZ(:,:,i),'FaceColor','k','FaceAlpha',0.7,'EdgeColor','k')
    surf(d3wakeX(:,:,i),d3wakeY(:,:,i),d3wakeZ(:,:,i),d3wakeG(:,:,i),'FaceAlpha',0.7)
end
plot3([0,0],[0,0],[0,dhubHeight],'k','LineWidth',2)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
zlim([0 inf])
c = colorbar('southoutside');
colormap jet
caxis([min(min(min(d3wakeG))) max(max(max(d3wakeG)))])
c.Label.String = 'Circulation \Gamma [m^2/s]';
view([-30 15])
daspect([1 1 1])

R = [];
end