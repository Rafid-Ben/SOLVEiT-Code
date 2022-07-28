function nBody_animate_new(Sdata, eventsData, nSteps, nShow, dim, p, view_spec,...
    title_opts, output, dt, nBodies, fraction, aniName, boxInfo)

boxInfoExists = true;
if nargin == 13
    boxInfoExists = false;
end

if view_spec == 'get'
[caz,cel] = view;
end

scale = 1e6; %scaling factor for the plot (currently microns); 1 for default
%scale = 1;


fac = 1.1; % amount to extend grid beyond largest value in each dimension
    dimx = fac*max(abs(Sdata(nSteps).pos(:,1)))*scale;
    dimy = fac*max(abs(Sdata(nSteps).pos(:,2)))*scale;
    dimz = fac*max(abs(Sdata(nSteps).pos(:,3)))*scale;
    dim_all = max([dimx dimy dimz]);
    minDimZ = min(Sdata(nSteps).pos(:,3))*scale * fac - (fac-1)*dimz;
if dim == 0
elseif dim == 1
    dimx = dim_all;
    dimy = dim_all;
    dimz = dim_all;
else
    dimx = dim*scale;
    dimy = dim*scale;
    dimz = dim*scale;
end
dimx = double(dimx);
dimy = double(dimy);
dimz = double(dimz);

close all
h = figure; %figure handle in order to be able to create a GIF file
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1.5;
ax.TickLabelInterpreter = 'latex';

axis equal
hold on

pos = (Sdata(1).pos) .* scale;

crashToggle = uint32(intmax('uint32'));
outOfBoundsToggle = uint32(0);

crashed_mon = (([eventsData(:).species] == 0) & ([eventsData(:).time] <= Sdata(1).time) & ...
    ([eventsData(:).event] == crashToggle));
outOfBounds_mon = (([eventsData(:).species] == 0) & ([eventsData(:).time] <= Sdata(1).time) & ...
    ([eventsData(:).event] == outOfBoundsToggle));

n_mon_crashed = sum(crashed_mon);
n_mon_outOfBounds = sum(outOfBounds_mon);
if (n_mon_crashed == 0)
    pos_crashed_mon = zeros(0, 3);
else
    pos_crashed_mon = zeros(n_mon_crashed, 3);
    j = 1;
    for i = 1:length(crashed_mon)
        if (crashed_mon(i) == 1)
            pos_crashed_mon(j,:) = eventsData(i).pos;
            j = j + 1;
        end
    end
end

if (n_mon_outOfBounds == 0)
    pos_outOfBounds_mon = zeros(0, 3);
else
    pos_outOfBounds_mon = zeros(n_mon_outOfBounds, 3);
    j = 1;
    for i = 1:length(outOfBounds_mon)
        if (outOfBounds_mon(i) == 1)
            pos_outOfBounds_mon(j,:) = eventsData(i).pos;
            j = j + 1;
        end
    end
end
pos_crashed_mon = pos_crashed_mon .* scale;
pos_outOfBounds_mon = pos_outOfBounds_mon .* scale;


mon = Sdata(1).species == 0;
pos_mon = (pos(mon,1:3));

dim = Sdata(1).species == 1;
pos_dim = (pos(dim,1:3));

tri = Sdata(1).species == 2;
pos_tri = (pos(tri,1:3));

neut = Sdata(1).species == 4;
pos_neut = (pos(neut,1:3));

monneg = Sdata(1).species == 5;
pos_monneg = (pos(monneg,1:3));

dimneg = Sdata(1).species == 6;
pos_dimneg = (pos(dimneg,1:3));

trineg = Sdata(1).species == 7;
pos_trineg = (pos(trineg,1:3));

n_mon = size(pos_mon,1);
n_dim = size(pos_dim,1);
n_tri = size(pos_tri,1);
n_neut = size(pos_neut,1);
n_monneg = size(pos_monneg,1);
n_dimneg = size(pos_dimneg,1);
n_trineg = size(pos_trineg,1);

%interval = 1; %only use when planning to display all steps

n_solv0 = n_dim + n_tri; % total number of solvated species injected

if (nSteps < 200)
   interval = 1; 
else
    interval = round(nSteps/200); % we only want to plot 100 distinct positions
end

color_crashed_mon = [0, 0, 1];
color_outOfBounds_mon = [1, 0, 1];

color_mon = [0.2 0.2 0.8]; %blue
color_dim = [0.9 0.2 0.2]; %red
color_tri = [0.2 0.9 0.2]; %green
color_neut = [0.9290, 0.6940, 0.1250];
color_monneg = [0 0 0];  %black
color_dimneg = [0.75, 0, 0.75]; %purple
color_trineg = [0, 0.75, 0.75]; %light blue


f_mon_crashed = plot3(pos_crashed_mon(1:n_mon_crashed,1),pos_crashed_mon(1:n_mon_crashed,2),pos_crashed_mon(1:n_mon_crashed,3),'o','Color',color_crashed_mon,'MarkerFaceColor',color_crashed_mon,'MarkerSize',12);
f_mon_outOfBounds = plot3(pos_outOfBounds_mon(1:n_mon_outOfBounds,1),pos_outOfBounds_mon(1:n_mon_outOfBounds,2),pos_outOfBounds_mon(1:n_mon_outOfBounds,3),'*','Color',color_outOfBounds_mon,'MarkerFaceColor',color_outOfBounds_mon,'MarkerSize',12);

f_mon = plot3(pos_mon(1:n_mon,1),pos_mon(1:n_mon,2),pos_mon(1:n_mon,3),'o','Color',color_mon,'MarkerFaceColor',color_mon,'MarkerSize',2);
f_dim = plot3(pos_dim(1:n_dim,1),pos_dim(1:n_dim,2),pos_dim(1:n_dim,3),'o','Color',color_dim,'MarkerFaceColor',color_dim,'MarkerSize',2);
f_tri = plot3(pos_tri(1:n_tri,1),pos_tri(1:n_tri,2),pos_tri(1:n_tri,3),'o','Color',color_tri,'MarkerFaceColor',color_tri,'MarkerSize',2);
f_neut = plot3(pos_neut(1:n_neut,1),pos_neut(1:n_neut,2),pos_neut(1:n_neut,3),'o','Color',color_neut,'MarkerFaceColor',color_neut,'MarkerSize',2);
f_monneg = plot3(pos_monneg(1:n_monneg,1),pos_monneg(1:n_monneg,2),pos_monneg(1:n_monneg,3),'o','Color',color_monneg,'MarkerFaceColor',color_monneg,'MarkerSize',2);
f_dimneg = plot3(pos_dimneg(1:n_dimneg,1),pos_dimneg(1:n_dimneg,2),pos_dimneg(1:n_dimneg,3),'o','Color',color_dimneg,'MarkerFaceColor',color_dimneg,'MarkerSize',2);
f_trineg = plot3(pos_trineg(1:n_trineg,1),pos_trineg(1:n_trineg,2),pos_trineg(1:n_trineg,3),'o','Color',color_trineg,'MarkerFaceColor',color_trineg,'MarkerSize',2);

if n_mon_crashed == 0
    f_mon_crashed = plot3(0, 0, 0,'o','Color',color_crashed_mon,'MarkerFaceColor',color_crashed_mon,'MarkerSize',1);
end
if n_mon_outOfBounds == 0
    f_mon_outOfBounds = plot3(0, 0, 0,'*','Color',color_outOfBounds_mon,'MarkerFaceColor',color_outOfBounds_mon,'MarkerSize',2);
end
if n_mon == 0
    f_mon = plot3(0,0,0,'o','Color',color_mon,'MarkerFaceColor',color_mon,'MarkerSize',2);
end
if n_dim == 0
    f_dim = plot3(0,0,0,'o','Color',color_dim,'MarkerFaceColor',color_dim,'MarkerSize',2);
end
if n_tri == 0
    f_tri = plot3(0,0,0,'o','Color',color_tri,'MarkerFaceColor',color_tri,'MarkerSize',2);
end
if n_neut == 0
    f_neut = plot3(0,0,0,'o','Color',color_neut,'MarkerFaceColor',color_neut,'MarkerSize',2);
end
if n_monneg == 0
    f_monneg = plot3(0,0,0,'o','Color',color_monneg,'MarkerFaceColor',color_monneg,'MarkerSize',2);
end
if n_dimneg == 0
    f_dimneg = plot3(0,0,0,'o','Color',color_dimneg,'MarkerFaceColor',color_dimneg,'MarkerSize',2);
end
if n_trineg == 0
    f_trineg = plot3(0,0,0,'o','Color',color_trineg,'MarkerFaceColor',color_trineg,'MarkerSize',2);
end

if (boxInfoExists)
    factor = length(boxInfo) / nSteps;
    [boxEdges, levelEdges] = plotBox(double(boxInfo(factor).minPoint .* scale),...
        double(boxInfo(factor).boxSize .* scale),double(boxInfo(factor).dims),double(boxInfo(factor).numLevel));
    boxPlot = plot3(boxEdges(:,1),boxEdges(:,2),boxEdges(:,3),'k');
    levelPlot = plot3(levelEdges(:,1),levelEdges(:,2),levelEdges(:,3),'b--');
end

if view_spec == 'reg'
elseif view_spec == 'x'
    view(90,0)
elseif view_spec == 'y'
    view(0,0)
elseif view_spec == 'z'
    view(0,90)
elseif view_spec == 'get'
    view(caz,cel)
end

%plottitle = title_opts;
%title(plottitle, 'FontSize', 20,'Interpreter','latex')
  
f.XDataSource = 'xi';
f.YDataSource = 'yi';
f.YDataSource = 'zi';
ax = gca;

set(ax, 'XLim', [-dimx dimx])
set(ax, 'YLim', [-dimy dimy])
set(ax, 'ZLim', [minDimZ dimz]) %2*dimz
hold on

xlabel ('x microns','Fontsize',20,'Interpreter','latex')
ylabel ('y microns','Fontsize',20,'Interpreter','latex')
zlabel ('z microns','Fontsize',20,'Interpreter','latex')

gif = aniName + '.gif'; %creating GIF file
video = VideoWriter(aniName + '.avi'); %creating video file
open (video);

xmin= -0.9*dimx; % left boundary of textbox
width = 0.74*dimx; % width of textbox
xtext = xmin + 0.02*dimx;

%plotcube2 (scale,xmin,width);
txt = text(xtext, 2.8*dimy, 1.9*dimz, num2str(dt));
%txt2 = text(xtext, 2.8*dimy, 1.8*dimz, num2str(nBodies/nSteps*1));
%txt3 = text(xtext, 2.8*dimy, 1.8*dimz, num2str(0));
%txt4 = text(xtext, 2.8*dimy, 1.7*dimz, num2str(0));
%txt5 = text(xtext, 2.8*dimy, 1.6*dimz, num2str(0));
%txt6 = text(xtext, 2.8*dimy, 1.5*dimz, num2str(0));
%txt7 = text(xtext, 2.8*dimy, 1.4*dimz, num2str(0));


frame = getframe(h); % setting up for .gif
[mov(:,:,1,1), map] = rgb2ind(frame.cdata, 256, 'nodither');

set(gcf,'color','w');


for i = 2:interval:nSteps
    
    crashed_mon = (([eventsData(:).species] == 0) & ([eventsData(:).time] <= Sdata(i).time) & ...
        ([eventsData(:).event] == crashToggle));
    outOfBounds_mon = (([eventsData(:).species] == 0) & ([eventsData(:).time] <= Sdata(i).time) & ...
        ([eventsData(:).event] == outOfBoundsToggle));

    n_mon_crashed = sum(crashed_mon);
    n_mon_outOfBounds = sum(outOfBounds_mon);
    if (n_mon_crashed == 0)
        pos_crashed_mon = zeros(0, 3);
    else
        pos_crashed_mon = zeros(n_mon_crashed, 3);
        j = 1;
        for k = 1:length(crashed_mon)
            if (crashed_mon(k) == 1)
                pos_crashed_mon(j,:) = eventsData(k).pos;
                j = j + 1;
            end
        end
    end

    if (n_mon_outOfBounds == 0)
        pos_outOfBounds_mon = zeros(0, 3);
    else
        pos_outOfBounds_mon = zeros(n_mon_outOfBounds, 3);
        j = 1;
        for k = 1:length(outOfBounds_mon)
            if (outOfBounds_mon(k) == 1)
                pos_outOfBounds_mon(j,:) = eventsData(k).pos;
                j = j + 1;
            end
        end
    end
    pos_crashed_mon = pos_crashed_mon .* scale;
    pos_outOfBounds_mon = pos_outOfBounds_mon .* scale;
    
    n_mon_crashed = size(pos_crashed_mon, 1);
    n_mon_outOfBounds = size(pos_outOfBounds_mon, 1);
    
    pos = (Sdata(i).pos) .* scale;
    
    mon = Sdata(i).species == 0; 
    pos_mon = (pos(mon,1:3));
    dim = Sdata(i).species == 1;
    pos_dim = (pos(dim,1:3));
    tri = Sdata(i).species == 2;
    pos_tri = (pos(tri,1:3));
    neut = Sdata(i).species == 4;
    pos_neut = (pos(neut,1:3));
    monneg = Sdata(i).species == 5; 
    pos_monneg = (pos(monneg,1:3));
    dimneg = Sdata(i).species == 6;
    pos_dimneg = (pos(dimneg,1:3));
    trineg = Sdata(i).species == 7;
    pos_trineg = (pos(trineg,1:3));
    
    n_mon = size(pos_mon,1);
    n_dim = size(pos_dim,1);
    n_tri = size(pos_tri,1);
    n_neut = size(pos_neut,1);
    n_monneg = size(pos_monneg,1);
    n_dimneg = size(pos_dimneg,1);
    n_trineg = size(pos_trineg,1);
    
    set(f_mon_crashed, 'XData', pos_crashed_mon(1:n_mon_crashed,1));
    set(f_mon_crashed, 'YData', pos_crashed_mon(1:n_mon_crashed,2));
    set(f_mon_crashed, 'ZData', pos_crashed_mon(1:n_mon_crashed,3));
    
    set(f_mon_outOfBounds, 'XData', pos_outOfBounds_mon(1:n_mon_outOfBounds,1));
    set(f_mon_outOfBounds, 'YData', pos_outOfBounds_mon(1:n_mon_outOfBounds,2));
    set(f_mon_outOfBounds, 'ZData', pos_outOfBounds_mon(1:n_mon_outOfBounds,3));
    
    set(f_mon, 'XData', pos_mon(1:n_mon,1));
    set(f_mon, 'YData', pos_mon(1:n_mon,2));
    set(f_mon, 'ZData', pos_mon(1:n_mon,3));
    
    
    set(f_dim, 'XData', pos_dim(1:n_dim,1));
    set(f_dim, 'YData', pos_dim(1:n_dim,2));
    set(f_dim, 'ZData', pos_dim(1:n_dim,3));
    
    set(f_tri, 'XData', pos_tri(1:n_tri,1));
    set(f_tri, 'YData', pos_tri(1:n_tri,2));
    set(f_tri, 'ZData', pos_tri(1:n_tri,3));
    
    set(f_neut, 'XData', pos_neut(1:n_neut,1));
    set(f_neut, 'YData', pos_neut(1:n_neut,2));
    set(f_neut, 'ZData', pos_neut(1:n_neut,3));
    
    
    set(f_monneg, 'XData', pos_monneg(1:n_monneg,1));
    set(f_monneg, 'YData', pos_monneg(1:n_monneg,2));
    set(f_monneg, 'ZData', pos_monneg(1:n_monneg,3));
    
    set(f_dimneg, 'XData', pos_dimneg(1:n_dimneg,1));
    set(f_dimneg, 'YData', pos_dimneg(1:n_dimneg,2));
    set(f_dimneg, 'ZData', pos_dimneg(1:n_dimneg,3));
    
    set(f_trineg, 'XData', pos_trineg(1:n_trineg,1));
    set(f_trineg, 'YData', pos_trineg(1:n_trineg,2));
    set(f_trineg, 'ZData', pos_trineg(1:n_trineg,3));
    
    if (boxInfoExists)
    %     [boxEdges, levelEdges] = plotBox(boxInfo(i*factor).minPoint .* scale,boxInfo(i*factor).boxSize .* scale,boxInfo(i*factor).dims,boxInfo(i*factor).numLevel);
        [boxEdges, levelEdges] = plotBox(double(boxInfo(i * factor).minPoint .* scale),...
            double(boxInfo(i * factor).boxSize .* scale),double(boxInfo(i * factor).dims),double(boxInfo(i * factor).numLevel));
    
        set(boxPlot, 'XData', boxEdges(:,1));
        set(boxPlot, 'YData', boxEdges(:,2));
        set(boxPlot, 'ZData', boxEdges(:,3));
        set(levelPlot,'XData',levelEdges(:,1));
        set(levelPlot,'YData',levelEdges(:,2));
        set(levelPlot,'ZData',levelEdges(:,3));
    end
    
    
%    nCrashed = size(find(Sdata(i).crash==1),1);
%    nFragment = size(find(Sdata(i).species==4),1);
    monomers = size(find(Sdata(end).species==0),1);
    dimers = size(find(Sdata(end).species==1),1);
    trimers = size(find(Sdata(end).species==2),1);
    monomersneg = size(find(Sdata(end).species==5),1);
    dimersneg = size(find(Sdata(end).species==6),1);
    trimersneg = size(find(Sdata(end).species==7),1);
%    percentFrag = nFragment/n_solv0*100;
    
    %set(txt,'String',['\textbf{Time: ' num2str(dt*1e9*i+dt*1e9, '%.0f') ' ns}'],'fontsize',12.5,'Interpreter','latex')
    set(txt,'String',['\textbf{Particles injected: ' num2str(round(nBodies/(nSteps*fraction))*i) '}'],'fontsize',12.5,'Interpreter','latex')
    %set(txt3,'String',['\textbf{Monomers: ' num2str(monomers) '}'],'fontsize',12.5,'Interpreter','latex')
    %set(txt4,'String',['\textbf{Dimers: ' num2str(dimers) '}'],'fontsize',12.5,'Interpreter','latex')
    %set(txt5,'String',['\textbf{Trimers: ' num2str(trimers) '}'],'fontsize',12.5,'Interpreter','latex')
    %legend('Mon','Dim','Tri','Neg Mon','Neg Dim','Neg Trim','fontsize',16,'Interpreter','latex','Location','northeast')
   

    drawnow update
    %pause(p);
    
    if output == 'G' %creating the GIF file
       
     %Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256);
      
      mov(:,:,1,i) = rgb2ind(frame.cdata, map, 'nodither');
    
        %writing GIF
        if i == 2 
          imwrite(imind,cm,gif,'gif', 'DelayTime', 0.01, 'Loopcount',inf); 
        else 
          imwrite(imind,cm,gif,'gif','WriteMode','append'); 
        end
        
    elseif output == 'V' %creating the video file
        
      %Capture the plot as an image 
      frame = getframe(h); 
      writeVideo(video,frame);
    
    else 
       fprintf ('Wrong input! Could not create video file or GIF.\n') 
       return
    end
        
end
% Create animated GIF
imwrite(mov, map, 'animation.gif', 'DelayTime', 0, 'LoopCount', inf);

close (video);

end

