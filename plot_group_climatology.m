%% read data
close all
types={'diurnal','diurnal_Om4','diurnal_Om025','dailyS0','dailyS0_Om4','dailyS0_Om025','diurnalspdx4','diurnalspdx025','dailyS0_T42','dailyS0_Om025_T42'};
Oms=[1,4,0.25,1,4,0.25,1,1,1,0.25];
sunmoves=[1,1,1,0,0,0,4,0.25,0,0];
lcs={'ko','ro','bo','k*','r*','b*','k^','kv','kd','bd'};

types={'diurnal','diurnal_Om025','dailyS0','dailyS0_Om025','diurnalspdx4','diurnalspdx025','dailyS0_T42','dailyS0_Om025_T42'};
Oms=[1,0.25,1,0.25,1,1,1,0.25];
sunmoves=[1,1,0,0,4,0.25,0,0];
lcs={'ko','ro','k*','r*','k^','kv','kd','rd'};

% types={'alb025','dailyS0'};
% Oms=[1,1];
% sunmoves=[0,0];
% lcs={'g*','k*'};

% types={'dailyS0','dailyS0_T42'};
% Oms=[1,1];
% sunmoves=[0,0];
% lcs={'k*','kd'};

S0=1360;
sigma=5.67e-8;
ntype=length(types);
plteachtype=false;
saveplt=false;

for itype=1:ntype
%for itype=3:3
cd(char(types(itype)))
lc=char(lcs(itype));
filehead='monthlyclimo_lonavg_';
Om=Oms(itype);
sunmove=sunmoves(itype);
tausws1=[0,0.5,1,2,4,6,8,12];
taulws1=[4];
[tausws,taulws]=meshgrid(tausws1,taulws1);
tausws=tausws(:); taulws=taulws(:);
chis=taulws./tausws;
ncase=length(tausws);


a=6371000;
g=9.81;
kappa=2/7;
P0=1000;
R=287.;
d2r=pi/180.;
f0=2*2*pi/86400;
for icase=1:ncase
    tausw=tausws(icase);
    taulw=taulws(icase);
    filename=[filehead,'s',erase(num2str(tausw),'.'),'l',num2str(taulw),'.nc'];
    if itype==1 && icase==1
        % read coordinate
        lat=ncread(filename,'lat'); nlat=length(lat);
        lev=ncread(filename,'lev'); nlev=length(lev);
        ilev=ncread(filename,'ilev'); nilev=length(ilev);
        iloglev=log(ilev);
        dp=(ilev(2:end)-ilev(1:end-1)).*100;  
        dlogp=iloglev(2:end)-iloglev(1:end-1);
        Th_T=(P0./lev').^kappa;
        
        % initialize matrices
        T_all=zeros(ncase,nlat,nlev);
        T_EQ_all=zeros(ntype,ncase,nlev);
        T_glob_all=zeros(ntype,ncase,nlev);
        Ts_EQ_all=zeros(ntype,ncase);
        U_all=zeros(ncase,nlat,nlev);
        U_EQ_all=zeros(ntype,ncase);
        U_EQ_tau1_all=zeros(ntype,ncase);
        U_EQ_mean_all=zeros(ntype,ncase);
        U_EQ_bot_all=zeros(ntype,ncase);
        V_all=zeros(ncase,nlat,nlev);
        VT_all=zeros(ncase,nlat,nlev);
        UV_all=zeros(ncase,nlat,nlev);
        OMEGAU_all=zeros(ncase,nlat,nlev);
        OMEGA_all=zeros(ncase,nlat,nlev);
        VT_intp_all=zeros(ncase,nlat);
        VT_intp_max_all=zeros(ntype,ncase);
        Psi_all=zeros(ncase,nlat,nlev);
        tdtrad_all=zeros(ncase,nlat,nlev);
        uv_all=zeros(ncase,nlat,nlev);
        vt_all=zeros(ncase,nlat,nlev);
        omegau_all=zeros(ncase,nlat,nlev);
        omegau_p_all=zeros(ncase,nlat,nlev);
        Psimax_all=zeros(ntype,ncase);
        Hadley_edge_all=zeros(ntype,ncase);
        S_all=zeros(ncase,nlat,nlev);
        S_glob_all=zeros(ncase,nlev);
        S_glob_tropmean_all=zeros(ntype,ncase);
        S_EQ_all=zeros(ntype,ncase);
        S_EQ_bot8_all=zeros(ntype,ncase);
        N2_EQ_all=zeros(ntype,ncase);
        U_EQ_index_all=zeros(ntype,ncase);
        omegau_glob_all=zeros(ntype,ncase);
        U_thwind_all=zeros(ncase,nlat,nlev);
        U_thwind_EQtop_all=zeros(ntype,ncase);
        Ptropp_EQ_all=zeros(ntype,ncase);
    end
    STRTind=[1,1,1]; % lat x lev x time
    CNTind=[Inf,Inf,Inf];
    STRDind=[1,1,1];
    if contains(char(types(itype)),'T42')
        STRTind=[2,1,1];
        STRDind=[2,1,1];
    end
    T_all(icase,:,:)=mean(ncread(filename,'T',STRTind,CNTind,STRDind),3);
    T_EQ_all(itype,icase,:)=mean(T_all(icase,nlat/2:nlat/2+1,:),2);
    T_glob_all(itype,icase,:)=mean(squeeze(T_all(icase,:,:)).*cosd(lat),1)/mean(cosd(lat));
    Ts_EQ_all(itype,icase)=mean(squeeze(T_all(icase,abs(lat)<3,nlev)));
    U_all(icase,:,:)=mean(ncread(filename,'U',STRTind,CNTind,STRDind),3);
    U_EQ_all(itype,icase)=(mean(squeeze(U_all(icase,abs(lat)<3,1))));
    U_EQ_tau1_all(itype,icase)=(mean(squeeze(U_all(icase,abs(lat)<3,[11:12])),[1,2]));
    U_EQ_mean_all(itype,icase)=mean(mean(squeeze(U_all(icase,abs(lat)<3,2:end))).*dlogp(2:end)')./mean(dlogp(2:end));
    U_EQ_bot_all(itype,icase)=mean(U_all(icase,abs(lat)<3,end-7));
    V_all(icase,:,:)=mean(ncread(filename,'V',STRTind,CNTind,STRDind),3);
    VT_all(icase,:,:)=mean(ncread(filename,'VT',STRTind,CNTind,STRDind),3);
    UV_all(icase,:,:)=mean(ncread(filename,'VU',STRTind,CNTind,STRDind),3);
    OMEGAU_all(icase,:,:)=mean(ncread(filename,'OMEGAU',STRTind,CNTind,STRDind),3);
    OMEGA_all(icase,:,:)=mean(ncread(filename,'OMEGA',STRTind,CNTind,STRDind),3);
    
    VT_intp_all(icase,:)=sum(squeeze(VT_all(icase,:,:)).*dp',2);
    V_intp_all(icase,:)=sum(squeeze(V_all(icase,:,:)).*dp',2);
    VT_intp_max_all(itype,icase)=max(abs(VT_intp_all(icase,:)));
    
    Psi_all(icase,:,:)=cumsum(squeeze(V_all(icase,:,:)).*dp',2).*(2*pi*a/g/1e9).*cosd(lat);
    [tmp,imax_lat]=max(Psi_all(icase,nlat/2+1:nlat/2+min(floor(nlat*0.4/Om),nlat/2),end-16:end));
    [Psimax_all(itype,icase),imax_lev]=max(tmp);
    imax_lat=imax_lat(imax_lev)+nlat/2;
    imax_lev=imax_lev+nlev-16-1;
    [tmp,imin_lat]=min(Psi_all(icase,nlat/2+1:end,1:end));
    [~,imin_lev]=min(tmp);
    imin_lat=imin_lat(imin_lev)+nlat/2;
    imin_lev=imin_lev+0;
    %if any(squeeze(Psi_all(icase,imax_lat:end,imax_lev))<0)
    %if any(squeeze(Psi_all(icase,nlat/2+1:end,imax_lev))<0)
    %if any(squeeze(Psi_all(icase,nlat/2+1:end,nlev-4))<0)
    %if any(squeeze(Psi_all(icase,nlat/2+1:end,nlev-10))<0)
    if any(squeeze(Psi_all(icase,nlat/2+1:end,imin_lev))<0)
        %Hadley_edge_all(itype,icase)=lat(min(find(squeeze(Psi_all(icase,imax_lat:end,imax_lev))<0))+(imax_lat-1));
        Hadley_edge_all(itype,icase)=lat(min(find(squeeze(Psi_all(icase,nlat/2+1:end,imax_lev))<Psimax_all(itype,icase)/20))+(nlat/2));
        %Hadley_edge_all(itype,icase)=lat(min(find(squeeze(Psi_all(icase,nlat/2+1:end,nlev-4))<0))+(nlat/2));
        %Hadley_edge_all(itype,icase)=lat(min(find(squeeze(Psi_all(icase,nlat/2+1:end,nlev-10))<0))+(nlat/2));
        %Hadley_edge_all(itype,icase)=lat(min(find(squeeze(Psi_all(icase,imax_lat:end,end-1))<0))+(imax_lat-1));
        %Hadley_edge_all(itype,icase)=lat(min(find(squeeze(Psi_all(icase,nlat/2+1:end,imin_lev))<0))+(nlat/2));
        Hadley_edge_all(itype,icase)=lat(min(find(squeeze(Psi_all(icase,nlat/2+1:end,imin_lev))<max(Psi_all(icase,nlat/2+1:end,imin_lev))/5))+(nlat/2));
        Hadley_edge_all(itype,icase)=lat(min(find(squeeze(Psi_all(icase,imax_lat:end,nlev-3))<Psimax_all(itype,icase)/4))+imax_lat-1);
        Hadley_edge_all(itype,icase)=lat(imax_lat);
        
        %[~,ilatmax]=max(Psi_all(icase,:,nlev-15));
        %Hadley_edge_all(itype,icase)=lat(ilatmax);
        
        Psiint=squeeze(sum(squeeze(Psi_all(icase,:,1:end)).*dp(1:end)',2)./sum(dp(1:end)));
        Psiint=mean(Psi_all(icase,:,1:end),3);
        ineg=min(find(Psiint(nlat/2+1:end)<max(Psiint)/20)); if length(ineg)==0; ineg=1;end
        Hadley_edge_all(itype,icase)=lat(ineg+nlat/2);
        
        %ilat_neg=max(min(find(any(Psi_all(icase,nlat/2+3:end,nlev-10:nlev)<0,3)))+nlat/2,nlat/2+1);
        %Hadley_edge_all(itype,icase)=lat(ilat_neg);
    else
        Hadley_edge_all(itype,icase)=NaN;
    end
    %Hadley_edge_all(itype,icase)=lat(imax_lat);
    tdtrad_all(icase,:,:)=mean(ncread(filename,'tdt_rad',STRTind,CNTind,STRDind),3).*86400;
    uv_all(icase,:,:)=UV_all(icase,:,:)-U_all(icase,:,:).*V_all(icase,:,:);
    vt_all(icase,:,:)=VT_all(icase,:,:)-T_all(icase,:,:).*V_all(icase,:,:);
    omegau_all(icase,:,:)=OMEGAU_all(icase,:,:)-OMEGA_all(icase,:,:).*U_all(icase,:,:);
    [omegau_p_all(icase,:,:),~]=gradient(squeeze(omegau_all(icase,:,:)),-lev.*100,lat);
    [Th_p,Th_y]=gradient(squeeze(T_all(icase,:,:)).*Th_T,lev.*100,lat.*d2r.*a);
    S_all(icase,:,:)=-Th_p./Th_T;
    S_glob_all(icase,:)=sum(squeeze(S_all(icase,:,:)).*cosd(lat))./sum(cosd(lat));
    S_glob_tropmean_all(itype,icase)=sum(S_glob_all(icase,:).*dp',2)./sum(dp);
    S_EQ_all(itype,icase)=S_all(icase,nlat/2,nlev);
    S_EQ_bot8_all(itype,icase)=S_all(icase,nlat/2,nlev-8);
    N2_EQ_all(itype,icase)=S_EQ_all(itype,icase)*P0*100*g^2/R/Ts_EQ_all(itype,icase)^2;
    U_EQ_index_all(itype,icase)=sum(squeeze(mean(omegau_p_all(icase,nlat/2:nlat/2+1,:),2)).*(dp./squeeze(mean(-OMEGA_all(icase,nlat/2:nlat/2+1,:),2))));
    U_thwind_all(icase,:,:)=cumsum(-(R/f0)./sind(lat).*(Th_y./Th_T).*(dp'./lev'./100),2,'reverse')./Om;
    U_thwind_EQtop_all(itype,icase)=max(mean(U_thwind_all(icase,nlat/2-5:nlat/2+6,:),2),[],3);
    Ptropp_EQ_all(itype,icase)=ilev(max(find((squeeze(S_all(icase,nlat/2+1,:)))>6e-4))+1);
end

if plteachtype
%% plot U
figure
logy=true;
polarmap(jet,0.7)
contours=[-240:20:240];
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(U_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    xticks([-90:30:90]); xlim([-90 90])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'U','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'U_levlat.png','png')
end

%% plot T
figure
logy=true;
polarmap(jet,0.7)
contours=[180:10:350];
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(T_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    xticks([-90:30:90]); xlim([-90 90])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'T','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'T_levlat.png','png')
end
%% plot V
% figure
% polarmap(jet,0.7)
% contours=[-1:0.1:1];
% for icase=1:ncase
%     subplot(length(taulws1),length(tausws1),icase)
%     contourf(lat,lev,squeeze(V_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
%     caxis([contours(1),contours(end)])
%     set(gca,'yscale','log','ydir','reverse','Fontsize',15)
%     xticks([-90:30:90]); xlim([-90 90])
%     yticks([0.06,0.3,1,3,10,30,100,300,900])
%     title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
%     if icase==1
%         ylabel('Pressure (mb)')
%         text(-190,0.02,'V','FontSize',20)
%     end
%     xlabel('latitude')
% end
% for ivert=1:length(taulws1)
%     hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
%     colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
% end
% set(gcf,'Position',[1         730        1680         218])

%% plot Psi
figure
polarmap(jet,0.7)
contours=[-10:1:10];
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(Psi_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    set(gca,'ydir','reverse','Fontsize',15)
    xticks([-90:45:90]); xlim([-90 90])
    %yticks([0.06,0.3,1,3,10,30,100,300,900])
    yticks([0:150:900])
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'\Psi','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'Psi_levlat.png','png')
end

%% plot Q_rad
figure
polarmap(jet,0.7)
contours=[-1:0.1:1];
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(tdtrad_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    set(gca,'ydir','reverse','Fontsize',15)
    xticks([-90:30:90]); xlim([-90 90])
    %yticks([0.06,0.3,1,3,10,30,100,300,900])
    yticks([0:150:900])
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'Q_{rad}','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'Rad_levlat.png','png')
end

%% plot VT
figure
polarmap(jet,0.7)
contours=[-100:10:100];
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(VT_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    set(gca,'ydir','reverse','Fontsize',15)
    xticks([-90:30:90]); xlim([-90 90])
    %yticks([0.06,0.3,1,3,10,30,100,300,900])
    yticks([0:150:900])
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'VT','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'VT_levlat.png','png')
end

%% plot v'T'
figure
polarmap(jet,0.7)
contours=[-20:2:20];
logy=false;
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(vt_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    xticks([-90:30:90]); xlim([-90 90])
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'vt','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'v_t_levlat.png','png')
end

%% plot u'v'
figure
polarmap(jet,0.7)
contours=[-50:5:50];
logy=true;
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(uv_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    xticks([-90:30:90]); xlim([-90 90])
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'uv','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'u_v_levlat.png','png')
end
%% plot u'omega'
figure
polarmap(jet,0.7)
contours=[-1:0.1:1]./100;
logy=true;
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(omegau_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    xticks([-90:30:90]); xlim([-90 90])
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'u\omega','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'u_omega_levlat.png','png')
end
%% plot u'omega'_p
figure
polarmap(jet,0.7)
contours=[-1:0.1:1]./2e5;
logy=true;
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(omegau_p_all(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    xticks([-90:30:90]); xlim([-90 90])
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'-u\omega_p','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'u_omega_p_levlat.png','png')
end

%% plot VT_intp
figure
cs=jet(ncase);
hold on
for icase=1:ncase
    plot(lat,squeeze(VT_intp_all(icase,:)),'Color',cs(icase,:),'LineWidth',2 ...
        ,'DisplayName',['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    xticks([-90:30:90]); xlim([-90 90])
    title('Northward Heat Transport')
    xlabel('latitude')
end
ss=plot(lat,lat.*0,'k--');
set(get(get(ss,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'FontSize',20)
set(gcf,'Position',[1   650   681   304])
legend('Location','NorthWestOutside')
if saveplt
saveas(gcf,'VTint_lat.png','png')
end
%% plot V_intp
figure
cs=jet(ncase);
hold on
for icase=1:ncase
    plot(lat,squeeze(V_intp_all(icase,:)),'Color',cs(icase,:),'LineWidth',2 ...
        ,'DisplayName',['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    xticks([-90:30:90]); xlim([-90 90])
    title('Northward Mass Transport')
    xlabel('latitude')
end
ss=plot(lat,lat.*0,'k--');
set(get(get(ss,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'FontSize',20)
set(gcf,'Position',[1   650   681   304])
legend('Location','NorthWestOutside')
if saveplt
saveas(gcf,'Vint_lat.png','png')
end

%% line plot
nplt=4;
figure
subplot(1,nplt,1)
scatter(tausws1,Psimax_all(itype,:),'ko')
xlabel('\tau_\infty^s')
title('\Psi_{max}')
set(gca,'FontSize',20)

subplot(1,nplt,2)
scatter(tausws1,VT_intp_max_all(itype,:),'ko')
xlabel('\tau_\infty^s')
title('Peak poleward heat trans')
set(gca,'FontSize',20)

subplot(1,nplt,3)
scatter(tausws1,U_EQ_all(itype,:),'ko')
xlabel('\tau_\infty^s')
title('Equatorial U max')
set(gca,'FontSize',20)

subplot(1,nplt,4)
scatter(tausws1,S_glob_tropmean_all(itype,:),'ko')
xlabel('\tau_\infty^s')
title('Static Stability')
set(gca,'FontSize',20)

set(gcf,'Position',[1         600        1350         350])
if saveplt
saveas(gcf,'Psi_VT_UEQ_S-s.png','png')
end

%% scatter plot
figure
nplt=3;

subplot(1,nplt,1)
scatter(S_glob_tropmean_all(itype,:),(Psimax_all(itype,:)))
xlabel('Static Stability')
ylabel('\Psi_{max}')
set(gca,'FontSize',20)

subplot(1,nplt,2)
scatter(S_glob_tropmean_all(itype,:),VT_intp_max_all(itype,:))
xlabel('Static Stability')
ylabel('Peak poleward heat trans')
set(gca,'FontSize',20)

subplot(1,nplt,3)
scatter(S_glob_tropmean_all(itype,:),U_EQ_all(itype,:))
xlabel('Static Stability')
ylabel('Equatorial U max')
set(gca,'FontSize',20)
set(gcf,'Position',[1         600        1350         350])
if saveplt
saveas(gcf,'Psi_VT_UEQ-S.png','png')
end

%% static stability
figure
logy=false;
cs=jet(ncase);
hold on
for icase=1:ncase
    plot(squeeze(S_glob_all(icase,:)),lev,'Color',cs(icase,:),'LineWidth',2 ...
        ,'DisplayName',['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    title('Static Stability')
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    ylabel('Pressure (mb)')
end
%set(gca,'xscale','log')
ylim([200 1000])
yticks([0:200:1000])
set(gca,'FontSize',20)
set(gcf,'Position',[1   650   681   304])
legend('Location','NorthWestOutside')
if saveplt
saveas(gcf,'Sglob_lev.png','png')
end
end
cd ..
end


%% radiative equilibrium temperature
lev1=10.^(0:0.01:3)';
nlev1=length(lev1);
Th_T1=(P0./lev1').^kappa;
Tskin=(S0*(cosd(lat)/pi)/(2*sigma)).^0.25;
Tskin0=(S0/pi/(2*sigma)).^0.25;
taulwcs_lev1=taulws./P0.*lev1';
taulwcs_lev=taulws./P0.*lev';
Trad1=zeros(ncase,nlat,nlev1);
Trad=zeros(ncase,nlat,nlev);
Tg=zeros(ncase,nlat);
noswabsorb=tausws==0;
n_noswabsorb=length(find(noswabsorb));
Trad1(noswabsorb,:,:)=reshape(Tskin,[1,nlat,1]).*reshape((1+taulwcs_lev1(noswabsorb,:)).^0.25,[n_noswabsorb,1,nlev1]);
Trad1(~noswabsorb,:,:)=reshape(Tskin,[1,nlat,1]).*reshape(((1+chis(~noswabsorb)+(1./chis(~noswabsorb)-chis(~noswabsorb)).*exp(-taulwcs_lev1(~noswabsorb,:)./chis(~noswabsorb)))).^0.25,[ncase-n_noswabsorb,1,nlev1]);
Trad(noswabsorb,:,:)=reshape(Tskin,[1,nlat,1]).*reshape((1+taulwcs_lev(noswabsorb,:)).^0.25,[n_noswabsorb,1,nlev]);
Trad(~noswabsorb,:,:)=reshape(Tskin,[1,nlat,1]).*reshape(((1+chis(~noswabsorb)+(1./chis(~noswabsorb)-chis(~noswabsorb)).*exp(-taulwcs_lev(~noswabsorb,:)./chis(~noswabsorb)))).^0.25,[ncase-n_noswabsorb,1,nlev]);
Tg(noswabsorb,:)=Tskin'.*(2+taulws(noswabsorb)).^0.25;
Tg(~noswabsorb,:)=Tskin'.*((1+chis(~noswabsorb))./2+(1-chis(~noswabsorb))./2.*exp(-taulws(~noswabsorb)./chis(~noswabsorb))).^0.25;

% solve for tropopause height: pr_tropp=p_tropp/ps
kappa=2/7;
pr_tropp=zeros(ncase,1);
for icase=1:ncase
    tau=taulws(icase);
    if tausws(icase)~=0
        chi=chis(icase);
        pr_tropp(icase)=fsolve(@(x) (1+chi-(chi-1/chi).*exp(-(tau/chi).*x)).*x.^(-4*kappa)-(1+chi+(1-chi)*exp(-tau/chi)),0.5);
    else
        pr_tropp(icase)=fsolve(@(x) (1+tau.*x).*x.^(-4*kappa)-(2+tau), 0.5);
    end
end
DP=max(1-pr_tropp,0).*P0;

% solve for thermal wind
Srad=Trad.*0; Urad_Om1=Trad.*0; Ugrad_EQ_Om1=squeeze(Urad_Om1(:,1,:)); Urad_Om1_taylor=Urad_Om1; DeltaH=zeros(ncase,1);DeltaV=zeros(ncase,1);lnth_tropp_to_surf=zeros(ncase,1);

for icase=1:ncase
% calculate thermal wind
Trad_y=-squeeze(Trad(icase,:,:))/4/a.*tand(lat);
[Thrad_p,Thrad_y]=gradient(squeeze(Trad(icase,:,:)).*Th_T,lev.*100,lat.*d2r.*a);
Srad(icase,:,:)=max(-Thrad_p./Th_T,0);
Urad_Om1(icase,:,:)=cumsum(-(R/f0)./sind(lat).*(Trad_y).*(dp'./lev'./100),2,'reverse');
Ugrad_EQ_Om1(icase,:)=fsolve(@(u) u+u.^2./(f0*a)./cosd(lat(nlat/2))-squeeze(Urad_Om1(icase,nlat/2,:)),squeeze(Urad_Om1(icase,nlat/2,:)));
seclat=1./cosd(lat); seclat(1)=0; seclat(end)=0;
tau_lev=taulws(icase).*lev./P0; dtau=taulws(icase)/P0.*dp./100;
Urad_Om1_taylor(icase,:,:)=(R*Tskin0/4/f0/a).*seclat.^(3/4).* ...
    (cumsum((1+tau_lev').^0.25./tau_lev'.*dtau',2,'reverse')...
    + cumsum((2-tau_lev'.^2)./tau_lev'./8./(1+tau_lev').^(3/4).*dtau',2,'reverse')./chis(icase));

% substract the wind below troppose height
idailymean1omega=find(strcmp(lcs,{'k*'}));
[~,iltropp]=min(abs(lev-Ptropp_EQ_all(idailymean1omega,icase)));
iltropp=iltropp-1;
ifric=min(abs(lev-710));
ibot=min(iltropp,ifric);
Urad_Om1(icase,:,:)=Urad_Om1(icase,:,:)-Urad_Om1(icase,:,ibot);
Urad_Om1(icase,:,ibot:end)=0;
Ugrad_EQ_Om1(icase,:)=Ugrad_EQ_Om1(icase,:)-Ugrad_EQ_Om1(icase,ibot);
Ugrad_EQ_Om1(icase,ibot:end)=0;
Urad_Om1_taylor(icase,:,:)=Urad_Om1_taylor(icase,:,:)-Urad_Om1_taylor(icase,:,ibot);
Urad_Om1_taylor(icase,:,ibot:end)=0;

% other useful quantities
DeltaH(icase)=std(Trad(icase,:,nlev));
[~,ilev1_tropp]=min(abs(lev1-pr_tropp(icase)*P0));
DeltaV(icase)=Trad1(icase,nlat/2,ilev1_tropp).*Th_T1(ilev1_tropp)-Trad1(icase,nlat/2,nlev1).*Th_T1(nlev1);
lnth_tropp_to_surf(icase)=log(Trad1(icase,nlat/2,ilev1_tropp).*Th_T1(ilev1_tropp))-log(Trad1(icase,nlat/2,nlev1).*Th_T1(nlev1));
end

%% critical tau for superrotation acceleration with taus
ftau=@(tau) (1/4).*log(((1+tau).^(1/4)-1)./((1+tau).^(1/4)+1))-atan((1+tau).^(1/4))-(1/10).*(1+tau).^(1/4).*(tau-4);
ftau4m=@(tau) ftau(4)-((1/4).*log(((1+tau).^(1/4)-1)./((1+tau).^(1/4)+1))-atan((1+tau).^(1/4))-(1/10).*(1+tau).^(1/4).*(tau-4));
taucrit4=fsolve(ftau4m,0.001)
pcrit4=taucrit4/4*1000
% plot ftau
figure
set(gca,'xscale','log')
hold on
fplot(ftau,[0 12],'k','LineWidth',2)
plot([0.01 4],[ftau(4) ftau(4)],'k--')
plot(4,ftau(4),'ko')
xlim([0.1,12])
ylim([-2.5,-1])
set(gca,'Fontsize',22)
xlabel('\tau')
xticks([1e-2,1e-1,1,2,4,8,12])

ylabel('$\mathcal{F}(\tau)$','Interpreter','latex')
set(gcf,'Position',[687   455   424   399])
if saveplt
saveas(gcf,'critical-tau.png','png')
end

%% plot Trad
if 1
figure
logy=true;
polarmap(jet,0.7)
contours=[180:10:350];
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(Trad(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    xticks([-90:30:90]); xlim([-90 90])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'T_{rad}','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'Trad_levlat.png','png')
end
end

%% plot Urad
if 1
figure
logy=true;
polarmap(jet,0.7)
contours=[-240:20:240];
for icase=1:ncase
    subplot(length(taulws1),length(tausws1),icase)
    contourf(lat,lev,squeeze(Urad_Om1(icase,:,:))',[-1000,contours,1000],'LineColor','none')
    caxis([contours(1),contours(end)])
    xticks([-90:30:90]); xlim([-90 90])
    if logy
        set(gca,'yscale','log','ydir','reverse','Fontsize',15)
        yticks([0.06,0.3,1,3,10,30,100,300,900])
    else
        set(gca,'ydir','reverse','Fontsize',15)
        yticks([0:150:900])
    end
    title(['\tau_\infty^s=',num2str(tausws(icase)),', \tau_\infty=',num2str(taulws(icase))])
    if icase==1
        ylabel('Pressure (mb)')
        text(-190,0.02,'U_{rad}','FontSize',20)
    end
    xlabel('latitude')
end
for ivert=1:length(taulws1)
    hp = get(subplot(length(taulws1),length(tausws1),(ivert*length(tausws1))-1),'Position');
    colorbar('Position', [hp(1)+hp(3)+0.12  hp(2)+0.02  0.01  hp(2)+hp(3)*7.7])
end
set(gcf,'Position',[1         730        1680         218])
if saveplt
saveas(gcf,'Urad_levlat.png','png')
end
end

%% scatter-all plot vs tausw
figure
nplt=4;
for itype=1:ntype
lc=char(lcs(itype));
subplot(1,nplt,2) % Held 2000
hold on
%yyaxis left
scatter(tausws1,max(Psimax_all(itype,:)-Psimax_all(itype,5),0).*Oms(itype).^2,100,lc)
xlabel('\tau_\infty^s')
title('\Psi_{max}')
set(gca,'FontSize',16)
set(gca,'yscale','log')
%set(gca,'xscale','log')
xticks([0 0.5 1 2 4 6 8 12])
if itype==ntype
    %fplot(@(x) 5.*x.^(-1.5),[0.5 10],'k-')
    %yyaxis right
    %plot(tausws1(1:4),DP(1:4).^1.5./150,'k--')
    %plot(tausws1(1:4),DP(1:4)./15,'k--')
    %plot(tausws1(1:4),DP(1:4).^1.5./S_all(1:4,16,23).^0.5./5000,'k--')
    %plot(tausws1(1:4),DP(1:4).^1.5./S_all(1:4,16,23).^0.5.*T_all(1:4,16,23).^3/2e11,'k--')
    %plot(tausws1(1:4),-log(pr_tropp(1:4)).*100,'k--')
    %plot(tausws1(1:4),(-log(pr_tropp(1:4))).^0.5./S_all(1:4,16,24).^0.5./2,'k--')
    %plot(tausws1(1:4),(-log(pr_tropp(1:4))).^1.*40,'k--')
    plot(tausws1(1:4),(-log(pr_tropp(1:4))).^1./Srad(1:4,16,24).^0.5./1.,'k--')
    plot(tausws1(1:4),(-log(pr_tropp(1:4))).^1./S_all(1:4,16,24).^0.5./1.,'k-')
    set(gca,'yscale','log')
    ylim([0.1 1e4])
    xlim([0 2])
end

subplot(1,nplt,1) % Held-Hou
hold on
%yyaxis left
scatter(tausws1,max(Psimax_all(itype,:)-Psimax_all(itype,5),0).*Oms(itype).^3,100,lc)
xlabel('\tau_\infty^s')
title('\Psi_{max}')
set(gca,'FontSize',16)
set(gca,'yscale','log')
%set(gca,'xscale','log')
xticks([0 0.5 1 2 4 6 8 12])
if itype==ntype
    %plot(tausws1(1:4),DP(1:4).^1.5./150,'k--')
    %plot(tausws1(1:4),DP(1:4)./15,'k--')
    %plot(tausws1(1:4),DP(1:4).^1.5./S_all(1:4,16,23).^0.5./5000,'k--')
    %plot(tausws1(1:4),DP(1:4).^1.5./S_all(1:4,16,23).^0.5.*T_all(1:4,16,23).^3/2e11,'k--')
    %plot(tausws1(1:4),-log(pr_tropp(1:4)).*100,'k--')
    %plot(tausws1(1:4),(-log(pr_tropp(1:4))).^2.5./Srad(1:4,16,24).^1./2,'k--')
    %plot(tausws1(1:4),(-log(pr_tropp(1:4))).^1.5.*DeltaH(1:4).^2.5./Srad(1:4,16,24).^1./2e5,'k--')
    plot(tausws1(1:4),(-log(pr_tropp(1:4))).^1.5.*DeltaH(1:4).^2.5./S_all(1:4,16,24).^1./2e5,'k-')
    plot(tausws1(1:4),(-log(pr_tropp(1:4))).^2.5.*DeltaH(1:4).^2.5./DeltaV(1:4)./2,'k--')
    set(gca,'yscale','log')
    ylim([0.1 1e4])
    xlim([0 2])
end

% subplot(1,nplt,2)
% hold on
% scatter(tausws1,Hadley_edge_all(itype,:).*Oms(itype).^0.5,100,lc)
% %set(gca,'yscale','log')
% if itype==ntype
%     %plot(tausws1(1:4),(-log(pr_tropp(1:4))).^0.33.*100,'k--')
%     %plot(tausws1(1:4),(-log(pr_tropp(1:4))).^0.5.*S_all(1:4,16,24).^0.5.*1.5e3,'k--')
%     plot(tausws1(1:4),(-log(pr_tropp(1:4))).^0.5.*DeltaV(1:4).^0.5.*1.5e3,'k--')
% end
% xlabel('\tau_\infty^s')
% xlim([0 2])
% title('Hadley cell edge')
% set(gca,'FontSize',16)

subplot(1,nplt,3)
hold on
showtop=0;showtau1=0;showmean=1;
if showtop
scatter(tausws1,U_EQ_all(itype,:).*Oms(itype),100,lc)
end
if showtau1
scatter(tausws1,U_EQ_tau1_all(itype,:).*Oms(itype),100,lc)
end
if showmean
scatter(tausws1,U_EQ_mean_all(itype,:).*Oms(itype),100,lc)
end
%plot(tausws1,U_thwind_EQtop_all(itype,:),[lc(1),'-'])
if itype==ntype
    if showtop
    plot(tausws1,Urad_Om1(:,nlat/2,1),'k--')
    plot(tausws1,Ugrad_EQ_Om1(:,1),'k-')
    %plot(tausws1,Urad_Om1_taylor(:,nlat/2,1),'k-.')
    end
    if showtau1
    plot(tausws1,squeeze(mean(Urad_Om1(:,nlat/2,11:12),3)),'k--')
    plot(tausws1,mean(Ugrad_EQ_Om1(:,11:12),2),'k-')
    plot(tausws1,squeeze(mean(Urad_Om1_taylor(:,nlat/2,11:12),3)),'k-.')
    end
    if showmean
    plot(tausws1,squeeze(mean(squeeze(Urad_Om1(:,nlat/2,2:end)).*dlogp(2:end)',2))./mean(dlogp(2:end)),'k--')
    plot(tausws1,squeeze(mean(Ugrad_EQ_Om1(:,2:end).*dlogp(2:end)',2))./mean(dlogp(2:end)),'k-')
    plot(tausws1,squeeze(mean(squeeze(Urad_Om1_taylor(:,nlat/2,2:end)).*dlogp(2:end)',2))./mean(dlogp(2:end)),'k-.')
    end
    
    xlabel('\tau_\infty^s')
    %set(gca,'yscale','log')
    xlim([0 12])
    %set(gca,'xscale','log')
    if showtop
    title('model top U_{eq}')
    end
    if showtau1
        title('U_{eq} at \tau=1')
    end
    if showmean
        title('vertically averaged U_{eq}')
    end
    set(gca,'FontSize',16)
    set(gcf,'Position',[1         600        1350         350])
end

subplot(1,nplt,4)
hold on
%scatter(tausws1,S_glob_all(itype,:),lc)
scatter(tausws1,S_EQ_all(itype,:),100,lc) % the one used
%scatter(tausws1,S_EQ_bot8_all(itype,:),100,lc)
if itype==ntype
    %plot(tausws1,mean(Srad(:,:,end),2),'k--')
    plot(tausws1,Srad(:,nlat/2,end),'k--') % the one used
    %plot(tausws1,Srad(:,nlat/2,end-8),'k--')
end
%scatter(tausws1,N2_EQ_all(itype,:),lc)
xlabel('\tau_\infty^s')
xlim([0 12])
title('Static Stability')
set(gca,'FontSize',16)
set(gcf,'Position',[1         600        1350         350])

if itype==ntype
%legend({'diurnal','diurnal,\Omega\times4','diurnal,\Omega/4','daily mean','daily mean,\Omega\times4','daily mean,\Omega/4','diurnal fast','diurnal slow','daily,high-res','daily,\Omega/4,high-res'},'Location','SouthEast')
legend({'diurnal','diurnal,\Omega/4','daily mean','daily mean,\Omega/4','diurnal fast','diurnal slow','daily,high-res','daily,\Omega/4,high-res'},'Location','SouthEast')
if saveplt
saveas(gcf,'./Psi_VT_UEQ_S-taus.png','png')
end
end
end

%% scatter-all Ptropp estimated vs diagnosed
figure
hold on
for itype=1:ntype
    lc=char(lcs(itype));
    scatter(Ptropp_EQ_all(itype,:),pr_tropp.*P0,200,lc)
end
xlabel('Diagnosed Tropopause Pressure (mb)')
ylabel('Estimated Tropopause Pressure (mb)')
set(gca,'FontSize',20)
xlim([400 1000])
ylim([400 1000])
plot([400 1000],[400 1000],'k--')
if saveplt
saveas(gcf,'Ptropp_diag-estim.png','png')
end

%% scatter-all Ptropp ~ exp(-tau)
figure
hold on
for itype=1:ntype
    lc=char(lcs(itype));
    scatter(exp(-tausws1),Ptropp_EQ_all(itype,:),100,lc)
    %if itype==7
    %text(exp(-tausws1*1.5),Ptropp_EQ_all(itype,:)-10,cellstr(strcat('\tau_\infty^s=',num2str(tausws1'))))
    %end
end
xlabel('exp(-\tau_\infty^s)')
%set(gca,'xscale','log')
%xticks(exp(-[10,5,0]))
%xticklabels(cellstr(strcat('\tau_\infty^s=',num2str([10,5,0]')))')
ylabel('Tropopause Pressure (mb)')
set(gca,'FontSize',16)
plot(exp(-tausws1),pr_tropp.*P0,'k--')
set(gcf,'Position',[411   614   300   317])
if saveplt
saveas(gcf,'pr_tropp-exptaus.png','png')
end

%% scatter-all log(P0/Ptropp) ~ tau
figure
hold on
for itype=1:ntype
    lc=char(lcs(itype));
    scatter(tausws1,log(P0./Ptropp_EQ_all(itype,:)),100,lc)
end
xlabel('\tau_\infty^s')
xlim([0 4])
ylabel('log(p_s/p_{tropp})')
set(gca,'FontSize',16)
plot(tausws1,log(1./pr_tropp),'k--')
set(gcf,'Position',[411   614   300   317])
if saveplt
saveas(gcf,'logpr_tropp-taus.png','png')
end


%% plot T-lev line
figure
hold on
ccs=jet(ncase);
set(gca,'yscale','log','ydir','reverse')
ylim([min(lev) max(lev)])
set(gca,'FontSize',18)
xlabel('Temperature (K)')
ylabel('Pressure (mb)')
for icase=1:ncase
    %plot(squeeze(Trad(icase,nlat/2,:)),lev,'Color',ccs(icase,:),'LineWidth',2)
    plot(mean(squeeze(Trad(icase,:,:)).*cosd(lat))/mean(cosd(lat)),lev,'Color',ccs(icase,:),'LineWidth',2)
end
for icase=1:ncase
for itype=1:ntype
    lc=char(lcs(itype));
    %scatter(squeeze(T_EQ_all(itype,icase,:)),lev,lc(2),'MarkerEdgeColor',ccs(icase,:))
    scatter(squeeze(T_glob_all(itype,icase,:)),lev,lc(2),'MarkerEdgeColor',ccs(icase,:))
end
end
%xlim([240 400])
xlim([250 400]-25)
legend(cellstr(strcat('\tau_\infty^s=',num2str(tausws1')))','Location','SouthEast')
set(gcf,'Position',[852   363   561   389])
if saveplt
saveas(gcf,'./T-lev.png','png')
end

%% plot Hadley boundary
figure
subplot(1,2,1) % Held Hou
hold on
for itype=1:ntype
lc=char(lcs(itype));
scatter(tausws1,Hadley_edge_all(itype,:).*Oms(itype).^(1),100,lc)
if itype==ntype
    plot(tausws1(1:4),(-log(pr_tropp(1:4))).^0.5.*30,'k--')
end
%set(gca,'yscale','log')
xlabel('\tau_\infty^s')
xlim([0 2])
ylim([0 50])
title('Held-Hou Hadley width (deg)')
set(gca,'FontSize',16)
end

subplot(1,2,2) %Held 2000
hold on
for itype=1:ntype
lc=char(lcs(itype));
scatter(tausws1,Hadley_edge_all(itype,:).*Oms(itype).^(1/3),100,lc)
if itype==ntype
    plot(tausws1(1:4),(lnth_tropp_to_surf(1:4)).^(1/3)*120,'k--')
end
xlabel('\tau_\infty^s')
xlim([0 2])
ylim([0 50])
title('Held-2000 Hadley width (deg)')
set(gca,'FontSize',16)
end
set(gcf,'Position',[731   639   829   316])
if saveplt
saveas(gcf,'./Hadley-boundary.png','png')
end


%% check resolution
% tausw_T42=0.5;taulw_T42=4;
% V_T42=squeeze(mean(ncread(['dailyS0_T42/monthlyclimo_lonavg_s',erase(num2str(tausw_T42),'.'),'l',erase(num2str(taulw_T42),'.'),'.nc'],'V'),3));
% lat_T42=ncread(['dailyS0_T42/monthlyclimo_lonavg_s',erase(num2str(tausw_T42),'.'),'l',erase(num2str(taulw_T42),'.'),'.nc'],'lat');
% Psi_T42=cumsum(squeeze(V_T42).*dp',2).*(2*pi*a/g/1e9).*cosd(lat_T42);
% 
% figure
% polarmap(jet,0.7)
% contours=[-10:1:10];
% 
% contourf(lat_T42,lev,Psi_T42',[-1000,contours,1000],'LineColor','none')
% caxis([contours(1),contours(end)])
% set(gca,'ydir','reverse','Fontsize',15)
% xticks([-90:30:90]); xlim([-90 90])
% %yticks([0.06,0.3,1,3,10,30,100,300,900])
% yticks([0:150:900])
% title(['\tau_\infty^s=',num2str(tausw_T42),', \tau_\infty=',num2str(taulw_T42)])
% if icase==1
%     ylabel('Pressure (mb)')
%     text(-190,0.02,'\Psi','FontSize',20)
% end
% xlabel('latitude')
% set(gcf,'Position',[374   311   408   436])
% set(gca,'fontsize',30)
% colorbar
% if saveplt
% saveas(gcf,'PsiT42_levlat.png','png')
% end


%% scatter-all plot vs static stability
% figure
% nplt=3;
% for itype=1:ntype
% lc=char(lcs(itype));
% subplot(1,nplt,1)
% hold on
% scatter(S_glob_tropmean_all(itype,:),Psimax_all(itype,:),lc)
% xlabel('Static Stability')
% title('\Psi_{max}')
% set(gca,'FontSize',20)
% if itype==ntype
%     legend({'diurnal','diurnal,\Omega\times4','diurnal,\Omega/4','daily mean','daily mean,\Omega\times4','daily mean,\Omega/4','diurnal fast','diurnal slow','daily mean,high-res'},'Location','NorthEast')
% end
% 
% subplot(1,nplt,2)
% hold on
% scatter(S_glob_tropmean_all(itype,:),Hadley_edge_all(itype,:),lc)
% xlabel('Static Stability')
% title('Hadley cell edge')
% set(gca,'FontSize',20)
% 
% subplot(1,nplt,3)
% hold on
% scatter(S_glob_tropmean_all(itype,:),U_EQ_all(itype,:),lc)
% xlabel('Static Stability')
% title('Equatorial U max')
% set(gca,'FontSize',20)
% set(gcf,'Position',[1         600        1350         350])
% 
% if itype==ntype
% if saveplt
% saveas(gcf,'./Psi_VT_UEQ-S.png','png')
% end
% end
% end