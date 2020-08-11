%% load files
load('/Users/manesh/Desktop/psychedelic_gradients/misc/fsa5_to_fsa10k/fsa_10k_midsurface_LR.mat')
mask=fsa_10k.mask;

% load('/Users/manesh/Desktop/psychedelic_gradients/comparative/DMT_PSILO_DMT_gradients_allPsychsMeanAligned.mat')
% load('/Users/manesh/Desktop/psychedelic_gradients/comparative/DMT/meanDMT_aligned/DMT_gradients_meanDMT_aligned.mat')
load('/Users/manesh/Desktop/psychedelic_gradients/comparative/margulies_template/DMT_PSILO_LSD_gradients_MargAligned.mat')


%% LSD behav

load('/Users/manesh/Desktop/psychedelic_gradients/FSA10k_results/ego_dissolution.mat')
load('/Users/manesh/Desktop/psychedelic_gradients/FSA10k_results/complex_imagery.mat')

ego=[ego_diss(:,1); ego_diss(:,2)];
img=[comp_imgry(:,1); comp_imgry(:,2)];

LSD_ego=ego_diss(:,1);
LSD_img=comp_imgry(:,1);

PLCB_ego=ego_diss(:,2);
PLCB_img=comp_imgry(:,2);

diff_ego=ego_diss(:,1) - ego_diss(:,2);
diff_img=comp_imgry(:,1) - comp_imgry(:,2);

%% Define terms for SurfStatLinMod

% Combined analysis

Group=[ones(17,1); zeros(17,1)];
Group = Group == 0;    

Group_fac = var2fac(Group, 'group'); % Group 0 = DMT, Group 1 = PLCB
Gr = term(Group_fac);

%% GLM looking for Group interactions between Cognitive scores and Gradient values
% 
% % ego dissolution 
% M                       = 1 +  Gr + E + E*Gr;
% 
% slm                     = SurfStatLinMod(Combined_G1', M, midsurf);
% contrast                = ego.*Gr.group0 - ego.*Gr.group1 ; %LSD > PLCB
% slm                     = SurfStatT(slm, contrast);
% 
% [pval peak clus clusid] = SurfStatP(slm,mask,0.05);
% ef                      = slm.t;
% ef(pval.C > 0.05)       = 0;
% 
% f                       = figure;
% SurfStatViewData        (ef, midsurf,'G Ego Diss Interaction (higher association in LSD at p<0.05')
% caxis                   ([ 0 max(slm.t) ])
% colormap                ([ 0.8 0.8 0.8 ; autumn  ])
% 
% % complex imagery
% 
% M                       = 1 +  Gr + I + I*Gr;
% 
% slm                     = SurfStatLinMod(Combined_G2', M, midsurf);
% contrast                = img.*Gr.group0 - img.*Gr.group1 ; %LSD > PLCB
% slm                     = SurfStatT(slm, -contrast);
% 
% [pval peak clus clusid] = SurfStatP(slm,mask,0.05);
% ef                      = slm.t;
% ef(pval.C > 0.05)       = 0;
% 
% f                       = figure;
% SurfStatViewData        (ef, midsurf,'G2 Complex Imagery Interaction; LSD>Placebo at p<0.05')
% caxis                   ([ 0 max(slm.t) ])
% colormap                ([ 0.8 0.8 0.8 ; autumn  ])


% Get values from clusters

% mean_clusterValue       = mean(Combined_G2(:,pval.C < 0.05),2);
% 
% f                       = figure; 
% SurfStatPlot            (mean_clusterValue, TotalCog, 1, Group_fac, 'LineWidth',2, 'MarkerSize',20);

%% LSD separately 

% Redefine terms (young)

% LSDE = term(LSD_ego);
% LSDI = term(LSD_img);

% LSDE = term(diff_ego);
% LSDI = term(diff_img);

E = term(ego);
I = term(img);
% IV = term(DMT_ImgVas');
% VImm = term(visImmersion');
% U= term(unity);


% GLM for ego dissolution

G2_diff=G1_DMT-G1_PCB;

M                       = 1 + E;
slm                     = SurfStatLinMod(G2_diff', M, fsa_10k);
slm                     = SurfStatT(slm, ego);

[pval peak clus clusid] = SurfStatP(slm,mask,0.05);
ef                      = slm.t;
ef(pval.C > 0.05)       = 0;

f                       = figure;
SurfStatViewData        (ef, fsa_10k, ' Ego Dissolution at p<0.05')
caxis                   ([ 0 max(slm.t) ])
colormap                ([ 0.8 0.8 0.8 ; autumn  ])

min(pval.C)


% GLM for complex imagery

G2_diff=G2_DMT-G2_PCB;

M                       = 1 + I;
slm                     = SurfStatLinMod(G2_diff', M, fsa_10k);
slm                     = SurfStatT(slm, img);

[pval peak clus clusid] = SurfStatP(slm,mask,0.05);
ef                      = slm.t;
ef(pval.C > 0.05)       = 0;

f                       = figure;
SurfStatViewData        (ef, fsa_10k, 'Complex Imagery at p<0.05')
caxis                   ([ 0 max(slm.t) ])
colormap                ([ 0.8 0.8 0.8 ; winter  ])


mean_clusterValue  = mean(G2_diff(pval.C < 0.05,:),1);

f=figure;
scatter                 (mean_clusterValue,ego, 85, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', [100/255, 188/255, 1])
xlabel('Mean Intra-Cluster Gradient Score')
ylabel('Complex Imagery DMT-Placebo')
title('G2 DMT-Placebo, Complex Imagery (p<0.05)')
ylim([0 1])
l1=lsline;
l1.Color='black';
l1.LineWidth=1;


