TNN-ADMM code for ``Denoising Electron Microscopy Data via a Hybrid Model Combining a Convolutional Neural Network and Tensor Nuclear Norm Minimization`` paper. WCVD and hybrid workflows are not included.

Run sample script: Exp_Script_Anscombe.m

- **Inputs:** `clean_crystal.mat` (ground truth) and `sm2_dec70.mat` (noisy video), included here. Each contains `data` in `[time, x, y]` order.
- **Parameters:** `lamin` controls regularization, `rhoin` the ADMM penalty, and `epsin` the stopping tolerance. `iterP.iM` sets the maximum iterations.
- **Outputs:** a comparison figure and a MAT file containing `Recon`, `rGTE_original`, `lamin`, and `rhoin`, saved in `Saved[Anscombe]`. 
