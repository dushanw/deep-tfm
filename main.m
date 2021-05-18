% 20200909 by Dushan N Wadduwage (wadduwage@fas.harvard.edu)
% 20201222 modified by DNW
% 20210518 modified by DNW

%% Notes and Requirements 
% (1).This code require CVX for optimization (http://cvxr.com/cvx/)
% (2).Beware that large-sized image reconstruction might be memory and time intensive.
% (3).To adapt to new data please update the sections "parameters" and "read and preprocess data"
% (4).This is a legacy version of the code and will is not maintained.
%     Please refer to https://github.com/dushanw/deep-tfm-vX.X 
%     (vX.X is the version, eg. the current version is v1.0) for the latest version.

%% system setup
clc; clear all; close all

% setup cvx
% run('../../_toolkits/cvx/cvx_startup.m')

%% parameters
pram.dx0          = 0.33;                 % original pixel size of image before resizeing to Nx Ny
pram.Nx           = 128;                  % 128x128 take 40mins per image on a computer with 1TB memory
pram.Ny           = 128;
pram.Nt           = 240;                  % #patterns
pram.gamma        = 5e-4;                 % weight of the wavelet sparsity prior in the optimization function
pram.wname        = 'db4';                % wavelet family used

pram.fileNameStem = sprintf('rec_Ny%d_Nx%d_Nt%d',pram.Ny,pram.Nx,pram.Nt);
pram.savepath     = ['./__results/' date '_deep-tfm_fig4/' pram.fileNameStem '/'];

mkdir(pram.savepath)                      % results will be saved to this folder

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
  Xhat_deep_noPr(:,:,:,i) = f_rec_inv_noPrior(pram,E,Y_deep(:,:,:,i),Y_tfm(:,:,:,i));         % no-prior 
  Xhat_deep_wlPr(:,:,:,i) = f_rec_inv_wlPrior(pram,E,Y_deep(:,:,:,i),pram.gamma,pram.wname);  % wavelet-prior with cvx
end
 save([pram.savepath 'reconstructed_' pram.fileNameStem '.mat'],'Xhat_deep_noPr','Xhat_deep_wlPr')
%load([pram.savepath 'reconstructed_' pram.fileNameStem '.mat'])

%% plot results
for i=1:size(Y_deep,4)
  figure('units','normalized','outerposition',[0 0 1 1])          
  imagesc([rescale(Y_tfm         (:,:,1,i)) ...
           rescale(Xhat_deep_noPr(:,:,1,i)) ...
           rescale(Xhat_deep_wlPr(:,:,1,i)) ]) ;axis image;colormap hot
  title('(1)Tfm-image (2)Deep-tfm-noPrior  (3)Deep-tfm-waveletPr')
  set(gca,'fontsize',24)

  saveas(gcf,[pram.savepath exp_names{i} '_fig.jpeg']);  
  close all
end
