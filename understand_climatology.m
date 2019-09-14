head='dailyS0_T42';
types={'','_damping','_damping_taupattern','_taupattern','_taupattern_albedo','_taupattern_albedo_humid','_obl23'};
names={'daily,high-res','damping','damping+\tau pattern','\tau pattern','\tau pattern+albedo','\tau pattern+albedo+moist','damping+obliquity'};

types={'','_taupattern','_taupattern_albedo','_taupattern_albedo_humid','_damping','_obl23'};
names={'daily,high-res','\tau pattern','\tau pattern+albedo','\tau pattern+albedo+moist','damping','damping+obliquity'};


ntype=length(types);
saveplt=true;

a=6371000;
g=9.81;
kappa=2/7;
P0=1000;
R=287.;
d2r=pi/180.;
f0=2*2*pi/86400;

for itype=1:ntype
filename=[head,char(types(itype)),'/monthlyclimo_lonavg_s0l4.nc'];
    if itype==1
        % read coordinate
        lat=ncread(filename,'lat'); nlat=length(lat);
        lev=ncread(filename,'lev'); nlev=length(lev);
        ilev=ncread(filename,'ilev'); nilev=length(ilev);
        iloglev=log(ilev);
        dp=(ilev(2:end)-ilev(1:end-1)).*100;  
        dlogp=iloglev(2:end)-iloglev(1:end-1);
        Th_T=(P0./lev').^kappa;
        
        Psi_all=zeros(ntype,nlat,nlev);
        U_all=zeros(ntype,nlat,nlev);
        T_all=zeros(ntype,nlat,nlev);
        Uthermal=zeros(ntype,nlat,nlev);
    end
    V=squeeze(mean(ncread(filename,'V'),3));
    U=squeeze(mean(ncread(filename,'U'),3));
    T=squeeze(mean(ncread(filename,'T'),3));
    Psi_all(itype,:,:)=cumsum(V.*dp',2).*(2*pi*a/g/1e9).*cosd(lat);
    U_all(itype,:,:)=U;
    T_all(itype,:,:)=T;
    [~,T_y]=gradient(squeeze(T),lev.*100,lat.*d2r.*a);
    Uthermal(itype,:,:)=cumsum(-(R/f0)./sind(lat).*(T_y).*(dp'./lev'./100),2,'reverse');
end

%% plot Psi
figure
polarmap(jet,0.7)
contours=[-10:1:10]*5;
for itype=1:ntype
    subplot(1,ntype,itype)
    contourf(lat,lev,squeeze(Psi_all(itype,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    set(gca,'ydir','reverse','Fontsize',15)
    xticks([-90:45:90]); xlim([-90 90])
    %yticks([0.06,0.3,1,3,10,30,100,300,900])
    yticks([0:150:900])
    title(char(names(itype)))
    if itype==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'\Psi','FontSize',20)
    end
    xlabel('latitude')
end

    hp = get(subplot(1,ntype,(ntype)-1),'Position');
    if ntype==7
    colorbar('Position', [hp(1)+hp(3)+0.13  hp(2)+0.1  0.01  hp(2)+hp(3)*6.5])
    end
    if ntype==6
    colorbar('Position', [hp(1)+hp(3)+0.15  hp(2)+0.1  0.01  hp(2)+hp(3)*5.2])
    end
    
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'understand_Psi_levlat.png','png')
end

%% plot U
figure
polarmap(jet,0.7)
contours=[-120:10:120];
logy=true;
for itype=1:ntype
    subplot(1,ntype,itype)
    contourf(lat,lev,squeeze(U_all(itype,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    set(gca,'ydir','reverse','Fontsize',15)
    xticks([-90:45:90]); xlim([-90 90])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    title(char(names(itype)))
    if itype==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'U','FontSize',20)
    end
    xlabel('latitude')
end

    hp = get(subplot(1,ntype,(ntype)-1),'Position');
    if ntype==7
    colorbar('Position', [hp(1)+hp(3)+0.13  hp(2)+0.1  0.01  hp(2)+hp(3)*6.5])
    end
    if ntype==6
    colorbar('Position', [hp(1)+hp(3)+0.15  hp(2)+0.1  0.01  hp(2)+hp(3)*5.2])
    end
    
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'understand_U_levlat.png','png')
end

%% plot Uthermal
figure
polarmap(jet,0.7)
contours=[-120:10:120];
logy=true;
for itype=1:ntype
    subplot(1,ntype,itype)
    contourf(lat,lev,squeeze(Uthermal(itype,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    set(gca,'ydir','reverse','Fontsize',15)
    xticks([-90:45:90]); xlim([-90 90])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    title(char(names(itype)))
    if itype==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'U_{th}','FontSize',20)
    end
    xlabel('latitude')
end

    hp = get(subplot(1,ntype,(ntype)-1),'Position');
    if ntype==7
    colorbar('Position', [hp(1)+hp(3)+0.13  hp(2)+0.1  0.01  hp(2)+hp(3)*6.5])
    end
    if ntype==6
    colorbar('Position', [hp(1)+hp(3)+0.15  hp(2)+0.1  0.01  hp(2)+hp(3)*5.2])
    end
    
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'understand_Uthermal_levlat.png','png')
end

%% plot T
figure
polarmap(jet,0.7)
contours=[200:2:320];
logy=true;
for itype=1:ntype
    subplot(1,ntype,itype)
    contourf(lat,lev,squeeze(T_all(itype,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    set(gca,'ydir','reverse','Fontsize',15)
    xticks([-90:45:90]); xlim([-90 90])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    title(char(names(itype)))
    if itype==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'T','FontSize',20)
    end
    xlabel('latitude')
end

    hp = get(subplot(1,ntype,(ntype)-1),'Position');
    if ntype==7
    colorbar('Position', [hp(1)+hp(3)+0.13  hp(2)+0.1  0.01  hp(2)+hp(3)*6.5])
    end
    if ntype==6
    colorbar('Position', [hp(1)+hp(3)+0.15  hp(2)+0.1  0.01  hp(2)+hp(3)*5.2])
    end
    
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'understand_T_levlat.png','png')
end