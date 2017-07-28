function m_plotMDS(plotcord,clabel)

%% visualize MDS generated point configuration

plotdim=size(plotcord,2);
% clabel=X(:,1);% marker color to show some property of objects
if nargin<2
    clabel=[0,0,1];
end
% set xlim
xratio=0.05; % space retained at two bounds
yratio=0.05;
xmin=min(plotcord(:,1));xmax=max(plotcord(:,1));
ymin=min(plotcord(:,2));ymax=max(plotcord(:,2));
xoffset=(xmax-xmin)*xratio;
yoffset=(ymax-ymin)*yratio;
xleft=xmin-xoffset;
xright=xmax+xoffset;
yleft=ymin-yoffset;
yright=ymax+yoffset;

hfig=figure('position',[100,100,800,500]);
Markersz=22;
if plotdim==2
    scatter(plotcord(:,1),plotcord(:,2),Markersz,clabel,'LineWidth',1,...
        'MarkerFaceColor','flat','MarkerEdgeColor',[0.75,0.75,0.75]);

elseif plotdim==3
    scatter3(plotcord(:,1),plotcord(:,2),plotcord(:,3),Markersz,clabel,...
        'LineWidth',1,'MarkerFaceColor','flat','MarkerEdgeColor',[0.75,0.75,0.75]);
    zratio=0.05;
    zmin=min(plotcord(:,3));zmax=max(plotcord(:,3));
    zoffset=(zmax-zmin)*zratio;
    zleft=zmin-zoffset;
    zright=zmax+zoffset;
    set(gca,'zlim',[zleft,zright],'zticklabel','');
end
%---**customize colormap**---
mcopy=hsv;% or 'jet', 'parula'
msub=mcopy(1:44,:);
mcolormap=flipud(msub);
colormap(hfig,mcolormap)
% colormap(hfig,jet)
% colorbar;

%-----classic code from myself--
set(gca,'Fontsize',13);
set(gca,'xlim',[xleft,xright],'xticklabel','');
set(gca,'ylim',[yleft,yright],'yticklabel','');

grid on 
box on
% set(gca,'GridLineStyle','--')
hold on
%----------------highlight some intersted locations
% load interestLoc % location and label 
% idim=size(interestLoc,2);
% if idim-1==2
%     x=interestLoc(:,1);
%     y=interestLoc(:,2);
%     ilabel=interestLoc(:,4);
%     sz=25;
%     scatter(x,y,sz,'Marker',s,'Markeredgecolor','k');
%     text(x,y,ilabel);
% elseif idim-1==3
%     x=interestLoc(:,1);
%     y=interestLoc(:,2);
%     z=interestLoc(:,3);
%     ilabel=interestLoc(:,4);
%     sz=25;
%     scatter3(x,y,z,sz,'Marker',s,'Markeredgecolor','k');
%     text(x,y,z,ilabel);
% end