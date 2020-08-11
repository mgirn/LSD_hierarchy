%% load directories
clear all
dir=genpath('surfstat');   
addpath(dir);
dir=genpath('BrainSpace-master');   
%% get files

parent_subfolders=deblank(icatb_listSubFolders('/lbc/lbc1/derivatives/meica_py27_v3.2_surface'));

filePattern_L1='*Rest1_FSA5_lh.mgh';
filePattern_R1='*Rest1_FSA5_rh.mgh';
filePattern_L2='*Rest3_FSA5_lh.mgh';
filePattern_R2='*Rest3_FSA5_rh.mgh';

LSD_sess1_dir = '/lbc/lbc1/derivatives/MG/psychedelic_data/LSD_functionals/Rest1/fsa5_mgh';
LSD_sess2_dir = '/lbc/lbc1/derivatives/MG/psychedelic_data/LSD_functionals/Rest2/fsa5_mgh';
LSD_data_L1 = deblank(icatb_listFiles_inDir(LSD_sess1_dir, filePattern_L1));
LSD_data_R1 = deblank(icatb_listFiles_inDir(LSD_sess1_dir, filePattern_R1));
LSD_data_L2 = deblank(icatb_listFiles_inDir(LSD_sess2_dir, filePattern_L2));
LSD_data_R2 = deblank(icatb_listFiles_inDir(LSD_sess2_dir, filePattern_R2));

PCB_sess1_dir = '/lbc/lbc1/derivatives/MG/psychedelic_data/PCB_functionals/Rest1/fsa5_mgh';
PCB_sess2_dir = '/lbc/lbc1/derivatives/MG/psychedelic_data/PCB_functionals/Rest2/fsa5_mgh';
PCB_data_L1 = deblank(icatb_listFiles_inDir(PCB_sess1_dir, filePattern_L1));
PCB_data_R1 = deblank(icatb_listFiles_inDir(PCB_sess1_dir, filePattern_R1));
PCB_data_L2 = deblank(icatb_listFiles_inDir(PCB_sess2_dir, filePattern_L2));
PCB_data_R2 = deblank(icatb_listFiles_inDir(PCB_sess2_dir, filePattern_R2));


%% compute LSD gradients
for k=1:15

    curr_file_L1=deblank(LSD_data_L1(k,:));
    curr_file_L2=deblank(LSD_data_L2(k,:));
    curr_file_R1=deblank(LSD_data_R1(k,:));
    curr_file_R2=deblank(LSD_data_R2(k,:));
    
    curr_Z = make_connectome_parcel_LSD(LSD_sess1_dir, LSD_sess2_dir, curr_file_L1, curr_file_R1, curr_file_L2, curr_file_R2);
    
    LSD_toRemove(k,:)=all(curr_Z==0,2); % get columns that contain a zero
   
    curr_remove=LSD_toRemove(k,:);

    curr_Z(curr_remove',:)=[];
    curr_Z(:,curr_remove)=[];
    
    gm = GradientMaps();
    
    curr_gm = gm.fit(curr_Z);
    
    gradients_LSD{k,:}.gradients = curr_gm.gradients;
    gradients_LSD{k,:}.lambda = curr_gm.lambda;
end
      
% save('/lbc/lbc1/derivatives/MG/psychedelic_gradients/FSA10k_results/gradients_LSD.mat','gradients_LSD');
% save('/lbc/lbc1/derivatives/MG/psychedelic_gradients/FSA10k_results/LSD_toRemove.mat','LSD_toRemove');
% 
% clear curr_Z LSD_toRemove gradients_LSD

%% compute placebo gradients

for k=1:15

    curr_file_L1=deblank(PCB_data_L1(k,:));
    curr_file_L2=deblank(PCB_data_L2(k,:));
    curr_file_R1=deblank(PCB_data_R1(k,:));
    curr_file_R2=deblank(PCB_data_R2(k,:));
    
    curr_Z = make_connectome_parcel_LSD(PCB_sess1_dir, PCB_sess2_dir, curr_file_L1, curr_file_R1, curr_file_L2, curr_file_R2);
    
    PCB_toRemove(k,:)=all(curr_Z==0,2);
    
    curr_remove=PCB_toRemove(k,:);
    
    curr_Z(curr_remove',:)=[];
    curr_Z(:,curr_remove)=[];
    
    gm = GradientMaps();
    
    curr_gm = gm.fit(curr_Z);
    
    gradients_PCB{k,:}.gradients = curr_gm.gradients;
    gradients_PCB{k,:}.lambda = curr_gm.lambda;
    
    disp(num2str(k))
    
    currgr=strcat('/lbc/lbc1/derivatives/MG/psychedelic_gradients/FSA10k_results/gradients_PCB_', num2str(k), '.mat');
    save(currgr,'gradients_PCB');
    
end

save('/lbc/lbc1/derivatives/MG/psychedelic_gradients/FSA10k_results/gradients_PCB_.mat','gradients_PCB');
save('/lbc/lbc1/derivatives/MG/psychedelic_gradients/FSA10k_results/PCB_toRemove_.mat','PCB_toRemove');

%% compute mean gradient



%% add back in missing values
cd /lbc/lbc1/derivatives/MG/psychedelic_gradients/FSA10k_results
oad('gradients_LSD_missingVals.mat')
load('gradients_PCB_missingVals.mat')
load('LSD_toRemove.mat')
load('PCB_toRemove.mat')

for i=1:15 % for each subject
%     LSD_gradient_vectors=gradients_LSD{i, 1}.gradients{:, 1};
    PCB_gradient_vectors=gradients_PCB{i, 1}.gradients{:, 1};
    
%      idx=find(LSD_toRemove(i,:));
    idx=find(PCB_toRemove(i,:));
    
    for k=1:10 % for each gradient
        
%         curr_gradient_LSD(:,k)=single(zeros(10000,1)); %vert x gradient x sub
        curr_gradient_PCB(:,k)=single(zeros(10000,1));
        
        c=0;
        for j=1:10000 % for each vertex
            
            if ~ismember(j,idx)
                c=c+1;
%                 curr_gradient_LSD(j,k)=LSD_gradient_vectors(c,k);
                curr_gradient_PCB(j,k)=PCB_gradient_vectors(c,k);
            else
%                 curr_gradient_LSD(j,k)=0;
                curr_gradient_PCB(j,k)=0;
            end
        end
    end
    
%     LSD_gradients_unaligned{1,i}=curr_gradient_LSD;
    PCB_gradients_unaligned{1,i}=curr_gradient_PCB;
end

% save('LSD_gradients_unaligned.mat','LSD_gradients_unaligned')
% save('PCB_gradients_unaligned.mat','PCB_gradients_unaligned')

%% align gradients

all_embeddings=[LSD_gradients_unaligned PCB_gradients_unaligned];

for j=1:30
    embeddings_array(:,:,j)=all_embeddings{1,j};
end
mean_embeddings=mean_LSDandPLCB.gradients{1,1};


[LSD_PCB_gradients.realigned, LSD_PCB_gradients.xfms] = mica_iterativeAlignment(all_embeddings,1,mean_embeddings);

save('LSD_PCB_gradients_aligned.mat','LSD_PCB_gradients')

%% display gradients

load('/lbc/lbc1/derivatives/MG/aging_gradients/fsa5_to_fsa10k/fsa_10k_midsurface_LR.mat')

SurfStatViewData(mean_PCB.gradients{1,1}(:,2),fsa_10k)
colormap(jet)
