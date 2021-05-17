% 20200909 by Dushan N Wadduwage (DNW)
% 20201222 modified by DNW
% test for real data

%% system setup
clc; clear all; close all

addpath(genpath('./_functionsAndLayers/'))
addpath('./_Datasets/')
addpath('./_ExtPatternsets/')

% setup cvx <cpy cvx here or give it as requirement for the code>
run('../_toolkits/cvx/cvx_startup.m')

%% parameters
pram.dx0     = 0.33;                % original pixel size of image before resizeing to Nx Ny
pram.Nx      = 16;
pram.Ny      = 16;
pram.Nt      = 16;
gamma        = 5e-4;
wname        = 'db4';
savepath     = ['./__results/' date '_tfm_mouse_20201224/'];
mkdir(savepath)

fileNameStem = sprintf('rec_Ny%d_Nx%d_Nt%d',pram.Ny,pram.Nx,pram.Nt);
disp(fileNameStem)

%% read and preprocess data
load('./Data_fig4.mat');

E = imresize(single(Data_fig4.patterns(:,:,1:pram.Nt)) ,[pram.Ny pram.Nx]);
E = E     -  mean(E    ,   3);
E = E     ./ max (E    ,[],3);

exp_names{1}    = 'surface';
exp_names{2}    = 'depth100um';
exp_names{3}    = 'depth200um';
exp_names{4}    = 'depth300um';
Y_deep(:,:,:,1) = imresize(single(Data_fig4.surface   (:,:,1:pram.Nt)),[pram.Ny pram.Nx]);
Y_deep(:,:,:,2) = imresize(single(Data_fig4.depth100um(:,:,1:pram.Nt)),[pram.Ny pram.Nx]);
Y_deep(:,:,:,3) = imresize(single(Data_fig4.depth200um(:,:,1:pram.Nt)),[pram.Ny pram.Nx]);
Y_deep(:,:,:,4) = imresize(single(Data_fig4.depth300um(:,:,1:pram.Nt)),[pram.Ny pram.Nx]);
Y_tfm (:,:,:,1) = imresize(single(Data_fig4.surface_wf   (:,:,1:end)) ,[pram.Ny pram.Nx]);
Y_tfm (:,:,:,2) = imresize(single(Data_fig4.depth100um_wf(:,:,1:end)) ,[pram.Ny pram.Nx]);
Y_tfm (:,:,:,3) = imresize(single(Data_fig4.depth200um_wf(:,:,1:end)) ,[pram.Ny pram.Nx]);
Y_tfm (:,:,:,4) = imresize(single(Data_fig4.depth300um_wf(:,:,1:end)) ,[pram.Ny pram.Nx]);

Y_deep  = Y_deep -  mean(Y_deep,3);
Y_deep  = Y_deep ./ max(max(Y_deep,[],1),[],2);

Y_tfm   = Y_tfm  ./ max(max(Y_tfm,[],1),[],2);

%% reconstruct and save results
for i=1:size(Y_deep,4)
  i
  Xhat_deep_noPr(:,:,:,i) = f_rec_inv_noPrior(pram,E,Y_deep(:,:,:,i),Y_tfm(:,:,:,i));  % no-prior 
  Xhat_deep_wlPr(:,:,:,i) = f_rec_inv_wlPrior(pram,E,Y_deep(:,:,:,i),gamma,wname);     % wavelet-prior with cvx
end
save([savepath 'reconstructed_' fileNameStem '.mat'],'Xhat_deep_noPr','Xhat_deep_wlPr')
%load([savepath 'reconstructed_' fileNameStem '.mat'])

%% plot results
for i=1:size(Y_deep,4)
  figure('units','normalized','outerposition',[0 0 1 1])          
  imagesc([rescale(Y_tfm         (:,:,1,i)) ...
           rescale(Xhat_deep_noPr(:,:,1,i)) ...
           rescale(Xhat_deep_wlPr(:,:,1,i)) ]) ;axis image;colormap hot
  title('(1)Tfm-image (2)Deep-tfm-noPrior  (3)Deep-tfm-waveletPr')
  set(gca,'fontsize',24)

  saveas(gcf,[savepath fileNameStem exp_names{i} '_fig.jpeg']);  
  close all
end
