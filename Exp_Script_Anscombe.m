clc
clear 
format long
format compact;
close all;

lamin = 30;  
rhoin = 0.5;
epsin = 1e-4; 
seed = 10;

%% Load data

True_struct = load('clean_crystal.mat');
True_raw = double(True_struct.data);

din = load('sm2_dec70.mat');
Noisy_raw = double(din.data);

%% Orientation

True_t  = permute(True_raw,  [2 3 1]);
Noisy_t = permute(Noisy_raw, [2 3 1]);


% Apply Anscombe to count data directly
A_noisy = anscombe(Noisy_t);

% Rescale transformed data for ADMM i.e., rescale to [0,1]
A_min = min(A_noisy(:));
A_max = max(A_noisy(:));
A_rng = A_max - A_min;
d = (A_noisy - A_min) / A_rng;

% Run TNN-ADMM on rescaled transformed data d
initD.T = d;
initD.A = @(X) X; 
initD.Y = d;


iterP = struct;
method = 'TNN-ADMM';
iterP.mode = 'unconstrained'; 
iterP.iM = 150;
iterP.outF = 5;
iterP.rho = rhoin;
iterP.lambda = lamin;
iterP.epsilon = epsin;  
[contD,AV] = tADMMv7(initD,iterP);

A_denoised = contD.X * A_rng + A_min;  % Undo ADMM scaling in Anscombe domain
lambda_hat_t = inv_anscombe(A_denoised); % Inverse Anscombe gives scaled up reconstruction

Recon = ipermute(lambda_hat_t, [2 3 1]);


%% Compute errors

rGTE_original =  tNormF(rescale(Recon) - rescale(True_raw)) / tNormF(rescale(True_raw))

%% Figures

fig = figure;
image_idx = 1;

subplot(1,3,1) 
imshow(squeeze(rescale(True_raw(image_idx,:,:)))', [0 1]); 
title(sprintf('Original; TNN = %.2e', TNNmap2(True_raw)));
colorbar;

subplot(1,3,2)
imshow(squeeze(rescale(Noisy_raw(image_idx,:,:)))', [0 1]); 
title(sprintf('With Poisson Noise; TNN = %.2e', TNNmap2(Noisy_raw)));
colorbar;

subplot(1,3,3)
imshow(squeeze(rescale(Recon(image_idx,:,:)))', [0 1]); 
title(sprintf('Denoised; TNN = %.2e; rGTE = %.2f',TNNmap2(Recon), rGTE_original));
colorbar;

sgtitle(sprintf('lambda = %g, rho = %g', lamin, rhoin));


outFolder = 'Saved[Anscombe]';
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end
fname2 = sprintf('_lambda_%g_rho_%g', lamin, rhoin);
filename = fullfile(outFolder, ['figure' fname2 '.fig']);
savefig(gcf, filename);

data_filename = fullfile(outFolder, ['data' fname2 '.mat']);
save(data_filename, 'Recon', 'rGTE_original', 'lamin', 'rhoin');


