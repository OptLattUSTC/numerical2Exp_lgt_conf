% ED code of 1D BHM for simulating the confinement-deconfinement transition
% under the mapping to the PXP model.
% for the case of a positive mass region
%
%
%  ---------------------------------
% defination in PXP model
%     01. atomic occupation vs gauge field, 
%         with lattice site index: k = 1,2,3,4,...,n 
%
%         index k |   ODD  |   EVEN
%         ------- | ------ | -------
%            ->   |    2   |   0/1
%            <-   |   0/1  |    2
%
%         where, -> = 0.5; <- = -0.5;
%         Notice: gauge field is fixed at the two edge sites {0,n+1}.
%
%     02. mapping to matter field, 
%         with index j = {[0,1],[1,2],...,[k-1,k],...,[n,n+1]}
%
%
%  ---------------------------------
%
%
%  ---------------------------------
%   Author: Wei-Yong Zhang
%   Email: zhangwya@mail.ustc.edu.cn
%          zhangwya@gmail.com
%   Revision: 1.0
%   Create Date: Sep. 24, 2022
%   Last Update: Nov. 11, 2022
%  ---------------------------------



%%
close all
clc

addpath(genpath('libs'))
addpath(genpath('includes'))


%% create a folder for saving the final results
save_path = [pwd,filesep,'analy_res_deconf'];
if ~exist(save_path,'dir')
    mkdir(save_path)
end


%% init state
L = 10;
st_init_up = zeros(1,L);
st_init_dn = [2 0 2 0 1 1 2 0 2 0];
% only for PXP mapping, original state without excitation
st_orig    = [2 0 2 0 2 0 2 0 2 0];

if (length(st_init_up) ~= L) || (length(st_init_dn) ~= L)
    error('Error! Invalid input parameters!')
end

%
N_up = sum(st_init_up);
N_dn = sum(st_init_dn);
% the maximum occupation number of up/dn state in each lattice site
nMax = 3; 


%% basis generate
basis = boson_basis_spinor(L,N_up,N_dn,nMax);
ns = basis.n_bs;
fprintf('Total basis number is %d.\n',ns)


%% update the basis to the PXP model
% state index search
basis =  map_bhm2PXP_state_Fcn(st_init_up,st_init_dn,basis,st_orig);
idx_init_st = basis.idx_init_st;

phi_init = zeros(ns,1);
phi_init(idx_init_st) = 1;


%% Hubbard parameters
J = 36.0 * 2 * pi;
U = 820.0 * 2 * pi;

w = 4*sqrt(2)*J^2/U;

%
U_uu = U;      % U_up_up
U_dd = U;      % U_dn_dn
U_ud = U;      % U_up_dn
J_up = J;
J_dn = J;

% staggered potential, 
% Notice: added to the odd lattice sites, if 'staG0>0'
staG0 = 0;  
Tilt_E = 5.7 * 220 * 2 * pi;
mu_up_lt = 0.0*2*pi*(1:L) + mod(1:L,2)*0;
mu_dn_lt = Tilt_E*(1:L) - mod(1:L,2)*staG0 * 2 * pi;

% boundary condition
BDC = 'obc';   % 'obc' - open boundary condition; 'pbs' - periodic case
% BDC = 'pbc';   


%% evolve time for quantum walk
nt = 201;
T = 98 * 1e-03;          % 100 ms
tl = linspace(0,T,nt);
dt = tl(2) - tl(1);


%% hamiltonian elements generate
tStart = tic; 
ham_elems = hamiltonian_1d_bhm2C_elements(basis,BDC);
tEnd = toc(tStart);
fprintf('Elapsed time of hamiltonian elements generation is %.6f seconds.\n',tEnd)


%% observables
% gauge_org_Gs = 0.5*ones(1,L);
% matter_org_Gs = zeros(1,L-1);
gauge_org_Gs = basis.gauge_st_pxp_orig;
matter_org_Gs = basis.matter_st_pxp_orig;


%% quench
s_q = map_bhm1d2PXP_quenchDynamics(basis,ham_elems,J_up,J_dn,...
    U_uu,U_dd,U_ud,mu_up_lt,mu_dn_lt,nt,tl,gauge_org_Gs,matter_org_Gs,...
    phi_init);


%% observables
psic_vs_tl = s_q.psic_vs_tl;
probl_vs_tl = s_q.probl_vs_tl;
density_up_Mt = s_q.density_up_Mt;
density_dn_Mt = s_q.density_dn_Mt;
new_obs_LGt = s_q.new_obs_LGt;
new_obs_LMt = s_q.new_obs_LMt;
density_Mt_GL = s_q.density_Mt_GL;
parity_Mt_GL = s_q.parity_Mt_GL;
G_corr_MtA = s_q.G_corr_MtA;
gauge_Field_Mt = s_q.gauge_Field_Mt;
gauge_abs_orderP_Mt = s_q.gauge_abs_orderP_Mt;
gauge_orderP_Mt = s_q.gauge_orderP_Mt;
GIPR = s_q.GIPR;
TransD = s_q.TransD;


%% save result
save([save_path,filesep,...
    sprintf('analy_deconf_delta@%02dHz.mat',staG0)],...
    'J','U','tl','staG0','Tilt_E','BDC',...
    'ham_elems','basis','T','nt',...
    'gauge_org_Gs',...
    'matter_org_Gs','psic_vs_tl','probl_vs_tl',...
    'density_up_Mt','density_dn_Mt','new_obs_LGt',...
    'new_obs_LMt','density_Mt_GL','parity_Mt_GL',...
    'G_corr_MtA','gauge_Field_Mt','gauge_abs_orderP_Mt',...
    'gauge_orderP_Mt','GIPR','TransD')


%%
figure('Color','w')
hold on
plot(tl*1000,new_obs_LGt,'-','LineWidth',2)
plot(tl*1000,new_obs_LMt,'-','LineWidth',2)
hold off
box on
ax = gca;
ax.FontSize = 14;
xlabel('Evolution time (ms)','FontSize',16)
% ylabel('\epsilon','FontSize',16)


%% density profile of the |up> state
% x = 1:L;
% y = tl*1000;
% 
% figure('Color','w','Position',[120 120 560 420])
% imagesc(x,y,density_up_Mt)
% colormap(hot)
% colorbar
% ax = gca;
% ax.FontSize = 14;
% xlabel('Sites','FontSize',16)
% ylabel('Evolution time (ms)','FontSize',16)


%% density profile of the |dn> state
x = 1:L;
y = tl*1000;

figure('Color','w','Position',[720 120 560 420])
% imagesc(x,y,density_dn_Mt)
imagesc(x,y,density_Mt_GL)
colormap(hot)
colorbar
ax = gca;
ax.FontSize = 14;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)


%%
figure('Color','w')
hold on
plot(tl*1000,gauge_abs_orderP_Mt,'-',...
    'LineWidth',2)
plot(tl*1000,gauge_orderP_Mt,'-',...
    'LineWidth',2)
hold off
box on
ax = gca;
ax.FontSize = 14;
xlabel('Evolution time (ms)','FontSize',16)
ylabel('\epsilon','FontSize',16)


%%
x = 1:L;
y = tl*1000;

% density profile of the |dn> state
figure('Color','w','Position',[720 120 520 300])
imagesc(x,y,density_Mt_GL)
colormap(hot)
colorbar
ax = gca;
ax.FontSize = 14;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)
% export_fig('-r300',sprintf('Density_profile_NumericSim_delta@%02dHz.png',staG0))


%
figure('Color','w','Position',[120 120 520 300])
imagesc(x(1:end)-1,y,G_corr_MtA)
% colormap(cMap_hot)
colormap(hot)
colorbar
ax = gca;
ax.FontSize = 14;
xlabel('\it r','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)
title('Correlation \itG(r,t)','FontSize',16)
% export_fig('-r300',sprintf('TwoPts_correlation_Deconf_NumericSim_delta@%02dHz.png',staG0))



%% for papers

x = 1:L;
y = tl*1000;

% density profile of the |dn> state
figure('Color','w','Position',[120 120 320 310])
imagesc(x,y,density_Mt_GL)
colormap(hot)
colorbar
ax = gca;
ax.FontSize = 16;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)
export_fig([save_path,filesep,...
    sprintf('Density_Numerical_delta@%02dHz.pdf',staG0)])


%
figure('Color','w','Position',[520 120 320 310])
imagesc(x(1:end-1)-1,y,G_corr_MtA(:,1:end-1))
% colormap(cMap_hot)
colormap(hot)
colorbar
ax = gca;
ax.FontSize = 16;
xlabel('\it r','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)
title('Correlation \itG(r,t)','FontSize',16)
export_fig([save_path,filesep,...
    sprintf('TwoPts_correlation_Numerical_delta@%02dHz.pdf',staG0)])


