%% Load files
load('/Users/manesh/Desktop/psychedelic_gradients/misc/fsa5_to_fsa10k/fsa_10k_midsurface_LR.mat')
mask=fsa_10k.mask;


%% LSD

LSD_gradients=DMT_PSILO_LSD_gradients_MargAligned.realigned(:,:,55:69);
LSD_PCB_gradients=DMT_PSILO_LSD_gradients_MargAligned.realigned(:,:,70:end);

G1_LSD = reshape(LSD_gradients(:,1,:),10000,15);
G2_LSD = reshape(LSD_gradients(:,2,:),10000,15);
G3_LSD = reshape(LSD_gradients(:,3,:),10000,15);

G1_PCB = reshape(LSD_PCB_gradients(:,1,:),10000,15);
G2_PCB = reshape(LSD_PCB_gradients(:,2,:),10000,15);
G3_PCB = reshape(LSD_PCB_gradients(:,3,:),10000,15);

Combined_G1 = [G1_LSD  G1_PCB];
Combined_G2 = [G2_LSD  G2_PCB];
Combined_G3 = [G3_LSD  G3_PCB];



%% Display gradient

f=figure; SurfStatViewData (mean(G1_LSD,2), fsa_10k, 'G1 LSD'); colormap (parula)
exportfigbo(f,'G1_LSD.png','png',10)

f=figure; SurfStatViewData (mean(G1_PCB,2), fsa_10k, 'G1 PCB'); colormap (parula)
exportfigbo(f,'G1_LSD_PCB.png','png',10)

f=figure; SurfStatViewData (mean(G2_LSD,2), fsa_10k, 'G2 LSD'); colormap (parula)
exportfigbo(f,'G2_LSD.png','png',10)

f=figure; SurfStatViewData (mean(G2_PCB,2), fsa_10k, 'G2 PCB'); colormap (parula)
exportfigbo(f,'G2_LSD_PCB.png','png',10)

f=figure; SurfStatViewData (mean(G3_LSD,2), fsa_10k, 'G3 LSD'); colormap (parula)
exportfigbo(f,'G3_LSD.png','png',10)

f=figure; SurfStatViewData (mean(G3_PCB,2), fsa_10k, 'G3 PCB'); colormap (parula)
exportfigbo(f,'G3_LSD_PCB.png','png',10)


%% Gradient Histograms

% 1st Gradient
G1_PCB_mean = mean(G1_PCB,2);
G1_LSD_mean = mean(G1_LSD,2);

f=figure; hold on;
P_h1                    = histogram(G1_PCB_mean(G1_PCB_mean~=0), 43, 'FaceColor', [0.55 0.55 0.55], 'EdgeAlpha', 0.5);
P_h2                    = histogram(G1_LSD_mean(G1_LSD_mean~=0), 43, 'FaceColor', [100/255, 188/255, 1], 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
title                   ('Principal Gradient (Blue = LSD)')
% xlim                    ([-0.18 0.1]);
% ylim                    ([0 4500]);
exportfigbo             (f,'G1_hist_LSD', 'png', 10)

% 2nd gradient
G2_PCB_mean = mean(G2_PCB,2);
G2_LSD_mean = mean(G2_LSD,2);

f=figure; hold on;
S_h1                    = histogram(G2_PCB_mean(G2_PCB_mean~=0), 43, 'FaceColor', [0.55 0.55 0.55], 'EdgeAlpha', 0.5);
S_h2                    = histogram(G2_LSD_mean(G2_LSD_mean~=0), 43, 'FaceColor', [100/255, 188/255, 1], 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
title                   ('Second Gradient (Blue = LSD)')
% xlim                    ([-0.13 0.13]);
% ylim                    ([0 3000]);
exportfigbo             (f,'G2_hist_LSD', 'png', 10)

% 3rd gradient
G3_PCB_mean = mean(G3_PCB,2);
G3_LSD_mean = mean(G3_LSD,2);

f=figure; hold on;
S_h1                    = histogram(G3_PCB_mean(G3_PCB_mean~=0), 43, 'FaceColor', [0.55 0.55 0.55], 'EdgeAlpha', 0.5);
S_h2                    = histogram(G3_LSD_mean(G3_LSD_mean~=0), 43, 'FaceColor', [100/255, 188/255, 1], 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
title                   ('Third Gradient (Blue = LSD)')
exportfigbo             (f,'G3_hist_LSD', 'png', 10)
% xlim                    ([-0.3 0.3]);
% ylim                    ([0 1900]);

%% Define surfstat terms to be used in the GLM
% Make group Logical
Group=[ones(15,1); zeros(15,1)];
Group = Group == 0;    %0 = DMT, 1 = PLCB

Group_fac = var2fac(Group, 'group'); % Group 0 = DMT, Group 1 = PLCB
Gr = term(Group_fac);

%% GLM group differences

%% Gradient 1
% PSILO > Placebo

M = 1 + Gr;

slm = SurfStatLinMod(Combined_G1', M, fsa_10k);
contrast                = Gr.group0 - Gr.group1;
slm                     = SurfStatT(slm, contrast);
[sig.pval sig.peak sig.clus sig.clusid] = SurfStatP(slm,mask,0.05);


ef_pos                      = slm.t;
ef_pos(sig.pval.C > 0.05)       = 0;

% f                       = figure;
% SurfStatViewData        (ef_pos, fsa_10k, 'LSD > Placebo')
% caxis                   ([ 0 max(slm.t) ])
% colormap                ([ 0.8 0.8 0.8 ; autumn ])

% sig.pval.thresh=0.05;
% 
% f=figure; SurfStatView(sig.pval, fsa_10k, 'RFT-FWE')


M = 1 + Gr;

slm                     = SurfStatLinMod(Combined_G1', M, fsa_10k);
contrast                = Gr.group1 - Gr.group0;
slm                     = SurfStatT(slm, contrast);
[sig.pval sig.peak sig.clus sig.clusid] = SurfStatP(slm,mask,0.05);

ef_neg                      = slm.t;
ef_neg(sig.pval.C > 0.05)       = 0;

% f                       = figure;
% SurfStatViewData        (ef_neg, fsa_10k, 'Old < Young G1')
% caxis                   ([ 0 max(slm.t) ])
% colormap                ([ 0.8 0.8 0.8 ; winter  ])

%% Display both

ef_both=ef_pos+(ef_neg*-1);
f                       = figure;
SurfStatViewData        (ef_both, fsa_10k, 'G1 LSD > Placebo (p<0.05, cdt<0.05)')
colormap                ([  flipud(winter); 0.8 0.8 0.8 ;autumn  ])
exportfigbo             (f,'G1_diff', 'png', 10)

%% Gradient 2
% LSD > Placebo

M = 1 + Gr;

slm                     = SurfStatLinMod(Combined_G2', M, fsa_10k);
contrast                = Gr.group0 - Gr.group1;
slm                     = SurfStatT(slm, contrast);
[sig.pval sig.peak sig.clus sig.clusid] = SurfStatP(slm,mask,0.01);

ef_pos                      = slm.t;
ef_pos(sig.pval.C > 0.05)       = 0;

% f                       = figure;
% SurfStatViewData        (ef_pos, fsa_10k, 'LSD > Placebo G2')
% caxis                   ([ 0 max(slm.t) ])
% colormap                ([ 0.8 0.8 0.8 ; autumn ])

% Old < Young

M = 1 + Gr;

slm                     = SurfStatLinMod(Combined_G2', M, fsa_10k);
contrast                = Gr.group1 - Gr.group0;
slm                     = SurfStatT(slm, contrast);
[sig.pval sig.peak sig.clus sig.clusid] = SurfStatP(slm,mask,0.005);

ef_neg                      = slm.t;
ef_neg(sig.pval.C > 0.05)       = 0;

% f                       = figure;
% SurfStatViewData        (ef_neg, fsa_10k, 'LSD < Placebo G2')
% caxis                   ([ 0 max(slm.t) ])
% colormap                ([ 0.8 0.8 0.8 ; winter  ])

%% Display both

ef_both=ef_pos+(ef_neg*-1);
f                       = figure;
SurfStatViewData        (ef_both, fsa_10k, 'G2 LSD > Placebo (p<0.05, cdt<0.005)')
colormap                ([  flipud(winter); 0.8 0.8 0.8 ;autumn  ])
exportfigbo             (f,'G2_diff', 'png', 10)

%% Gradient 3

% Old < Young

M                       = 1 + Gr;

slm                     = SurfStatLinMod(Combined_G3', M, fsa_10k);
contrast                = Gr.group0 - Gr.group1;
slm                     = SurfStatT(slm, contrast);
[pval peak clus clusid] = SurfStatP(slm,mask,0.005);

ef_pos                     = slm.t;
ef_pos(pval.C > 0.05)       = 0;

% f                       = figure;
% SurfStatViewData        (ef, fsa_10k, 'LSD > Placebo G3')
% caxis                   ([ 0 max(slm.t) ])
% colormap                ([ 0.8 0.8 0.8 ; winter ])
% % exportfigbo             (f,'/host/yeatman/local_raid/alex/Manesh_aging/Results/Old<Young_G3', 'png', 10)

% Old > Young

M                       = 1 + Gr;

slm                     = SurfStatLinMod(Combined_G3', M, fsa_10k);
contrast                = Gr.group1 - Gr.group0;
slm                     = SurfStatT(slm, contrast);
[pval peak clus clusid] = SurfStatP(slm,mask,0.005);

ef_neg                      = slm.t;
ef_neg(pval.C > 0.05)       = 0;

% f                       = figure;
% SurfStatViewData        (ef, fsa_10k, 'LSD > Placebo G3')
% caxis                   ([ 0 max(slm.t) ])
% colormap                ([ 0.8 0.8 0.8 ; autumn  ])
% % exportfigbo             (f,'/host/yeatman/local_raid/alex/Manesh_aging/Results/Old>Young_G3', 'png', 10)



%% Display both

ef_both=ef_pos+(ef_neg*-1);
f                       = figure;
SurfStatViewData        (ef_both, fsa_10k, 'G3 LSD > Placebo (p<0.05, cdt<0.005)')
colormap                ([  flipud(winter); 0.8 0.8 0.8 ;autumn  ])
exportfigbo             (f,'G3_diff', 'png', 10)

%% Parcellate into Yeo Networks and display differences in mean gradient values on Spider plot

% cd                      /host/yeatman/local_raid/alex/Manesh_aging

% yeo                     = load('/Users/manesh/Desktop/Gradients/yeo_fsaverage5.mat');
% yeo                     = yeo.yeo;
load('/Users/manesh/Desktop/psychedelic_gradients/FSA10k_results/Archive/yeo10k.mat');
yeo = yeo10k';

DANidx                                  = yeo==2; 
FPNidx                                  = yeo==3;
DMNidx                                  = yeo==4;
VISidx                                  = yeo==5;
LIMidx                                  = yeo==6;
SMidx                                   = yeo==7;
SALidx                                  = yeo==8; 
net_toRemove=[1 6];
network.names           = {'DAN','FPCN','DN','VIS','SM', 'SAL'};

% load(Users/manesh/Desktop/psychedelic_gradients/FSA10k_results/Archive/yeo10k.mat');
% load('/Users/manesh/Desktop/aging_gradients/LATEST_FSA10k/gradients_10k/Figures/network_labels/yeo10k_DNsep.mat')


% network.names           = {'DAN','FPCN','VIS','LIM','SM','SAL','DNmtl', 'DNcore', 'DNsub3'};
% net_toRemove=[1 4 9];

% Gradient 1

[yFull_G1 ~]            = mica_surfData2parcelData(Combined_G1', yeo10k);
yFull_G1(:,net_toRemove)=[];

% get t vals and df for significance

M = 1 + Gr;

slm                     = SurfStatLinMod(yFull_G1, M);
contrast                = Gr.group0 - Gr.group1;  % DMT > Placebo
slm                     = SurfStatT(slm, contrast);

%%

DMT                    = find(Gr.group0==1);
PLCB                   = find(Gr.group1==1);

networks_PLCB_G1       = yFull_G1(PLCB,:);
networks_DMT_G1         = yFull_G1(DMT,:);

m_networks_PLCB_G1      = mean(networks_PLCB_G1,1);
m_networks_DMT_G1        = mean(networks_DMT_G1,1);
std_networks_PLCB_G1    = std(networks_PLCB_G1,1);

networks_PLCB_zeros = [0 0 0 0 0 0];

for i=1:size(networks_DMT_G1,1)    %every sub
    for k=1:size(networks_DMT_G1,2)    %every network
        normalized_DMT_G1(i,k)=(networks_DMT_G1(i,k)-m_networks_PLCB_G1(k))/std_networks_PLCB_G1(k);
    end
end
m_normalized_DMT_G1=mean(normalized_DMT_G1,1);


% m_Young_G1              = mean(yFull_G1(young,:),1);

f                       = figure;
spider_test             ([m_normalized_DMT_G1',networks_PLCB_zeros'], [], repmat([-2.5 1.5],6,1), ...
    network.names,[.8 0 0; 0 0 0],gca ) ;
set(gcf,'color','w');
exportfigbo             (f,'G1_spider_DMT', 'png', 10)  % Old = black, Red = young

% T STATS               = -0.9113   -2.7774   -3.0729   -1.3642    2.2876   -1.7163    1.1607    3.6757

%% Gradient 2

[yFull_G2 ~]            = mica_surfData2parcelData(Combined_G2', yeo10k);
yFull_G2(:,net_toRemove)=[];
% 
M                       = 1 + Gr;

slm                     = SurfStatLinMod(yFull_G2, M);
contrast                = Gr.group0 - Gr.group1;
slm                     = SurfStatT(slm, contrast);

%% 
DMT                     = find(Gr.group0==1);
PLCB                   = find(Gr.group1==1);

networks_PLCB_G2       = yFull_G2(PLCB,:);
networks_DMT_G2         = yFull_G2(DMT,:);
m_networks_PLCB_G2       = mean(yFull_G2(PLCB,:),1);
m_networks_DMT_G2         = mean(yFull_G2(DMT,:),1);
std_networks_PLCB_G2    = std(yFull_G2(PLCB,:),1);
networks_PLCB_zeros = [0 0 0 0 0 0];

for i=1:size(networks_DMT_G2,1)
    for k=1:size(networks_DMT_G2,2)
        normalized_DMT_G2(i,k)=(networks_DMT_G2(i,k)-m_networks_PLCB_G2(k))/std_networks_PLCB_G2(k);
    end
end
m_normalized_DMT_G2=mean(normalized_DMT_G2,1);


spider_test             ([m_normalized_DMT_G2',networks_PLCB_zeros'], [], repmat([-2.5 1.5],6,1), ...
    network.names,[.8 0 0; 0 0 0],gca ) ; %[-0.09 0.08]
set(gcf,'color','w');
exportfigbo             (f,'G2_spider_DMT', 'png', 10) 

% T STATS               = -1.7889    1.1378    2.2704    1.2213    1.7362    2.6353   -1.1799   -3.3834

%% Gradient 3

[yFull_G3 ~]            = mica_surfData2parcelData(Combined_G3', yeo10k);
yFull_G3(:,net_toRemove)=[];
% 
M                       = 1 + Gr;

slm                     = SurfStatLinMod(yFull_G3, M);
contrast                = Gr.group0 - Gr.group1;
slm                     = SurfStatT(slm, contrast);

%% 
DMT                     = find(Gr.group0==1);
PLCB                   = find(Gr.group1==1);

networks_PLCB_G3       = yFull_G3(PLCB,:);
networks_DMT_G3         = yFull_G3(DMT,:);
m_networks_PLCB_G3       = mean(yFull_G3(PLCB,:),1);
m_networks_DMT_G3         = mean(yFull_G3(DMT,:),1);
std_networks_PLCB_G3    = std(yFull_G3(PLCB,:),1);
networks_PLCB_zeros = [0 0 0 0 0 0];

for i=1:size(networks_DMT_G3,1)
    for k=1:size(networks_DMT_G3,2)
        normalized_DMT_G3(i,k)=(networks_DMT_G3(i,k)-m_networks_PLCB_G3(k))/std_networks_PLCB_G3(k);
    end
end
m_normalized_DMT_G3=mean(normalized_DMT_G3,1);


spider_test             ([m_normalized_DMT_G3',networks_PLCB_zeros'], [], repmat([-2.5 1.5],6,1), ...
    network.names,[.8 0 0; 0 0 0],gca ) ; %[-0.09 0.08]
set(gcf,'color','w');
exportfigbo             (f,'G3_spider_DMT', 'png', 10) 

%% Generate csv files for R joy plot

% load yeo annotation file

[~, lh_annot_lab_yeo, lh_annot_ctb_yeo] = fs_read_annotation('/Users/manesh/Desktop/Gradients/LATEST/Figures/OLD/nonCT_controlled/joy/to_make/lh.Yeo2011_7Networks_N1000.annot');
[~, rh_annot_lab_yeo, rh_annot_ctb_yeo] = fs_read_annotation('/Users/manesh/Desktop/Gradients/LATEST/Figures/OLD/nonCT_controlled/joy/to_make/rh.Yeo2011_7Networks_N1000.annot');
load('/Users/manesh/Desktop/psychedelic_gradients/FSA10k_results/Archive/yeo_annotations/lh_annot_lab_yeo10k.mat')
load('/Users/manesh/Desktop/psychedelic_gradients/FSA10k_results/Archive/yeo_annotations/rh_annot_lab_yeo10k.mat')

yeoL                                    = zeros(1,length(lh_annot_lab_yeo));
yeoR                                    = zeros(1,length(rh_annot_lab_yeo));

netLabels                               = lh_annot_ctb_yeo.table(:,5);

for i                                  = 1:length(netLabels)
    indicesL                            = lh_annot_lab_yeo == netLabels(i);
    yeoL(indicesL)                      = i;
    
    indicesR                            = rh_annot_lab_yeo == netLabels(i);
    yeoR(indicesR)                      = i;
end

yeo                                     = [yeoL,yeoR];



% Find vertices which contribute to each network

VISidx                                  = yeo==2;
SMidx                                   = yeo==3;
DANidx                                  = yeo==4;
SALidx                                  = yeo==5;
LIMidx                                  = yeo==6;
FPidx                                   = yeo==7;
DMNidx                                  = yeo==8;

%% Produce csv files for import into R for joy plots

% Mean the gradients across subjects

mean_DMT_G1 = mean(G1_DMT,2);
mean_PCB_G1 =mean(G1_PCB,2);

mean_DMT_G2 = mean(G2_DMT,2);
mean_PCB_G2 = mean(G2_PCB,2);


% Create the csv for R

% cd /host/yeatman/local_raid/alex/Manesh_aging/Results

% Principle Gradient

DMT_principle_yeo(:,1)                = mean_DMT_G1;    %10000 x 1
DMT_principle_yeo(:,2)                = yeo';
ONEidx                                  = DMT_principle_yeo(:,2)==1;
DMT_principle_yeo(ONEidx,:)           = [];
csvwrite                                ('DMT_G1_yeo.csv',DMT_principle_yeo);


PLCB_principle_yeo(:,1)                  = mean_PCB_G1;
PLCB_principle_yeo(:,2)                  = yeo';
ONEidx                                  = PLCB_principle_yeo(:,2)==1;
PLCB_principle_yeo(ONEidx,:)             = [];
csvwrite                                ('PLCB_G1_yeo.csv',PLCB_principle_yeo);

% Second Gradient

G2_DMT_yeo(:,1)                   = mean_DMT_G2;
G2_DMT_yeo(:,2)                   = yeo';
ONEidx                                  = G2_DMT_yeo(:,2)==1;
G2_DMT_yeo(ONEidx,:)              = [];
csvwrite                                ('DMT_G2_yeo.csv',G2_DMT_yeo);


G2_PLCB_yeo(:,1)                     = mean_PCB_G2;
G2_PLCB_yeo(:,2)                     = yeo';
ONEidx                                  = G2_PLCB_yeo(:,2)==1;
G2_PLCB_yeo(ONEidx,:)                = [];
csvwrite                                ('PLCB_G2_yeo.csv',G2_PLCB_yeo);

% Third Gradient

G3_young_yeo(:,1)                    = G3_PCB_mean';
G3_young_yeo(:,2)                    = yeo';
ONEidx                                  = G3_young_yeo(:,2)==1;
G3_young_yeo(ONEidx,:)               = []
csvwrite                                ('Young_G3_yeo.csv',G3_young_yeo);


G3_old_yeo(:,1)                      = G3_DMT_mean';
G3_old_yeo(:,2)                      = yeo';
ONEidx                                  = G3_old_yeo(:,2)==1;
G3_old_yeo(ONEidx,:)                 = []
csvwrite                                ('Old_G3_yeo.csv',G3_old_yeo);

%% Plot gradients on scatter (Manifold)

blue=[70 130 180]/256; %SM
purple=[120 18 134]/256; %VIS
green=[0 118 14]/256; %DAN
magenta=[196 58 250]/256; %SAL
cream=[220 248 164]/256; %LIM
orange=[230 148 34]/256; %FPN
red=[205 62 78]/256; % DN


yellow=[255 255 0]/256;
darkblue=[0 0 130]/256;
maroon=[205 62 78]/256;
grey=[0.8 0.8 0.8];

% 
% load('/Users/manesh/Desktop/psychedelic_gradients/FSA10k_results/Archive/yeo10k.mat');
% yeo = yeo10k';
% VISidx                                  = yeo==5;
% SMidx                                   = yeo==7;
% DANidx                                  = yeo==2; %8 
% SALidx                                  = yeo==8; %2
% LIMidx                                  = yeo==6;
% FPNidx                                  = yeo==3;
% DMNidx                                  = yeo==4;

load('/Users/manesh/Desktop/aging_gradients/LATEST_FSA10k/gradients_10k/Figures/network_labels/yeo10k_DNsep.mat')
yeo = yeo10k';
VISidx                                  = yeo==5;
SMidx                                   = yeo==7;
DANidx                                  = yeo==2;
SALidx                                  = yeo==8;
LIMidx                                  = yeo==6;
FPNidx                                  = yeo==3;

DNcoreidx                               = yeo==16;
DNmtlidx                                = yeo==15;
DNsub3idx                               = yeo==17;

% 
for i=1:10000
    colors{i}=grey;
end

for i=find(VISidx)
    colors{i}=purple;
end

for i=find(SMidx)
    colors{i}=blue;
end

for i=find(DANidx)
    colors{i}=green;
end

for i=find(SALidx)
    colors{i}=magenta;
end

% for i=find(LIMidx)
%     colors{i}=cream;
% end

for i=find(FPNidx)
    colors{i}=orange;
end

for i=find(DNcoreidx)
    colors{i}=red; %yellow
end

for i=find(DNmtlidx)
    colors{i}=red; %darkblue
end

for i=find(DNsub3idx)
    colors{i}=red; %maroon
end
% 

% for i=find(DNidx)
%     colors{i}=red;
% end


% 
% G2_diff=mean(G2_DMT,2)-mean(G2_PLCB,2);
% 
% 
% for i=find(G2_diff<0)'
%  colors{i}=blue;
% end
% 
% for i=find(G2_diff>0)'
%  colors{i}=red;
% end



colors_array=[];
for i=1:10000
    colors_array=vertcat(colors_array, colors{1,i});
end

f=figure
for i=1:9
    s=scatter(G2_LSD(:,i), G1_LSD(:,i),6,colors_array, 'filled');
    ylim([-0.16 0.16])
    xlim([-0.16 0.16])
    hold on
end
xlabel('G2', 'FontSize', 12);
ylabel('G1', 'FontSize', 12);
title('Gradient G1-G2 Manifold LSD');
exportfigbo(f,'manifold_LSD', 'png', 10) 

f=figure
for i=1:9
    
    s=scatter(G2_PCB(:,i), G1_PCB(:,i),6,colors_array, 'filled');
    ylim([-0.16 0.16])
    xlim([-0.16 0.16])
    hold on
end
xlabel('G2', 'FontSize', 12);
ylabel('G1', 'FontSize', 12);
title('Gradient 1-2 Manifold LSD Placebo');
exportfigbo(f,'manifold_LSD_PCB', 'png', 10) 



%%

% G1_neg10=sort(Combined_G1(:,i),'ascend');
% G1_neg10=mean(G1_neg10(1:ceil(length(G1_neg10)*0.1))); 

G1_diff=G1_DMT-G1_PCB;

for i=1:17
%     G1diff_NP(i,1)=abs(median(Combined_G1(:,i)-min(Combined_G1(:,i))));
%     G1diff_NP(i,2)=abs(median(Combined_G1(:,i)-max(Combined_G1(:,i))));
%     G1_diffDMTPLCB(i,2)=abs(max(Combined_G1(:,i+15)-min(Combined_G1(:,i+15))));
    G1_range_DMTPCB(i,1)=min(G1_DMT(:,i)-max(G1_PCB(:,i)));
end

for i=1:30
    G2diff_NP(i,1)=abs(median(Combined_G2(:,i)-min(Combined_G2(:,i))));
    G2diff_NP(i,2)=abs(median(Combined_G2(:,i)-max(Combined_G2(:,i))));
end



DMT_PCB_G1_NEG=[G1diff_NP(16:30,1) G1diff_NP(1:15,1)];
DMT_PCB_G1_POS=[G1diff_NP(16:30,2) G1diff_NP(1:15,2)];
DMT_PCB_G1_NP=[DMT_PCB_G1_NEG DMT_PCB_G1_POS];
save('DMT_PCB_G1_NP.mat','DMT_PCB_G1_NP')

DMT_PCB_G2_NEG=[G2diff_NP(16:30,1) G2diff_NP(1:15,1)];
DMT_PCB_G2_POS=[G2diff_NP(16:30,2) G2diff_NP(1:15,2)];
DMT_PCB_G2_NP=[DMT_PCB_G2_NEG DMT_PCB_G2_POS];
save('DMT_PCB_G2_NP.mat','DMT_PCB_G2_NP')


[~,P_G1(1),~,STATS] = ttest2(G1diff_NP(16:30,1), G1diff_NP(1:15,1)); % PCB-DMT G1 negative span
[~,P_G1(2),~,STATS] = ttest2(G1diff_NP(16:30,2), G1diff_NP(1:15,2)); % PCB-DMT G1 positive span

[~,P_G2(1),~,STATS] = ttest2(G2diff_NP(16:30,1), G2diff_NP(1:15,1)); % PCB-DMT G1 negative span
[~,P_G2(2),~,STATS] = ttest2(G2diff_NP(16:30,2), G2diff_NP(1:15,2)); % PCB-DMT G1 positive span

[~,P_G1,~,STATS] = ttest(G1_diffLP(:,1), G1_diffLP(:,2)); % PCB-DMT G1 positive span


%% 

G1_mean=mean([G1_DMT, G1_PCB],2);
G1_med=median(G1_mean);

G1_mean2=G1_mean;
G1_mean2(G1_mean<G1_med)=0;
G1_unimodal_idx=find(G1_mean2);

G1_mean3=G1_mean;
G1_mean3(G1_mean>G1_med)=0;
G1_transmodal_idx=find(G1_mean3);

G1_diffMEAN=mean(G1_diff,2);


boxplot([G1_diffMEAN(G1_unimodal_idx), G1_diffMEAN(G1_transmodal_idx)])

scatter(G1_diffMEAN(G1_unimodal_idx))













