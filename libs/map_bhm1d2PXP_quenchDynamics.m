function [s] = map_bhm1d2PXP_quenchDynamics(basis,ham_elems,J_up,J_dn,...
    U_uu,U_dd,U_ud,mu_up_lt,mu_dn_lt,nt,tl,gauge_org_Gs,matter_org_Gs,...
    phi_init)


%%
s = struct();


%%
ham_cur = hamiltonian_1d_bhm2C(basis,ham_elems,J_up,J_dn,...
    U_uu,U_dd,U_ud,mu_up_lt,mu_dn_lt);

ham = ham_cur.ham;
fprintf('\nNumber of nonzero hamiltonian elements: %d.\n',nnz(ham))


%% observables
psic_vs_tl = [];
probl_vs_tl = [];

density_up_Mt = [];
density_dn_Mt = [];

new_obs_LGt = [];
new_obs_LMt = [];

density_Mt_GL = [];
parity_Mt_GL = [];
G_corr_MtA = [];
gauge_Field_Mt = [];

gauge_abs_orderP_Mt = [];
gauge_orderP_Mt = [];

GIPR = [];
TransD = [];


%%
dt = tl(2)-tl(1);
psic = phi_init;
probl = abs(conj(psic).*psic);

stat_nC = map_bhm1d2PXP_statistics_Fcn(psic,basis,gauge_org_Gs,matter_org_Gs);
density_up_Mt = cat(1,density_up_Mt,stat_nC.density_up_ltC);
density_dn_Mt = cat(1,density_dn_Mt,stat_nC.density_dn_ltC);
density_Mt_GL = cat(1,density_Mt_GL,stat_nC.density_ltC_GL);
parity_Mt_GL = cat(1,parity_Mt_GL,stat_nC.parity_ltC_GL);
G_corr_MtA = cat(1,G_corr_MtA,stat_nC.G_dis_r);
gauge_Field_Mt = cat(1,gauge_Field_Mt,stat_nC.gauge_DensitylC);
new_obs_LGt = cat(1,new_obs_LGt,stat_nC.observablEs);
new_obs_LMt = cat(1,new_obs_LMt,stat_nC.observablMs);

gauge_abs_orderP_Mt = cat(1,gauge_abs_orderP_Mt,...
    stat_nC.gauge_abs_orderP_Cur);
gauge_orderP_Mt = cat(1,gauge_orderP_Mt,...
    stat_nC.gauge_orderP_Cur);

A = stat_nC.G_dis_r;
GIPR = cat(1,GIPR,sum(A.^2) / sum(A)^2);
TransD = cat(1,TransD,sum(A .* ((1:length(A))-1)));

psic_vs_tl = cat(2,psic_vs_tl,psic);
probl_vs_tl = cat(2,probl_vs_tl,probl);

%
tStart = tic; 
for kk = 1:nt-1
    fprintf('Current process: %04d / %04d.\n',kk,nt)
    
    tSC = tic;
    psic = expv(-1i*dt,ham,psic,1e-7,30);
    tEC = toc(tSC);
    fprintf('time for evolution: %.6f seconds.\n',tEC)
    
    probl = abs(conj(psic).*psic);
      
    stat_nC = map_bhm1d2PXP_statistics_Fcn(psic,basis,gauge_org_Gs,matter_org_Gs);
    density_up_Mt = cat(1,density_up_Mt,stat_nC.density_up_ltC);
    density_dn_Mt = cat(1,density_dn_Mt,stat_nC.density_dn_ltC);
    density_Mt_GL = cat(1,density_Mt_GL,stat_nC.density_ltC_GL);
    parity_Mt_GL = cat(1,parity_Mt_GL,stat_nC.parity_ltC_GL);
    G_corr_MtA = cat(1,G_corr_MtA,stat_nC.G_dis_r);
    gauge_Field_Mt = cat(1,gauge_Field_Mt,stat_nC.gauge_DensitylC);
    new_obs_LGt = cat(1,new_obs_LGt,stat_nC.observablEs);
    new_obs_LMt = cat(1,new_obs_LMt,stat_nC.observablMs);
    
    gauge_abs_orderP_Mt = cat(1,gauge_abs_orderP_Mt,...
        stat_nC.gauge_abs_orderP_Cur);
    gauge_orderP_Mt = cat(1,gauge_orderP_Mt,...
        stat_nC.gauge_orderP_Cur);
    
    A = stat_nC.G_dis_r;
    GIPR = cat(1,GIPR,sum(A.^2) / sum(A)^2);
    TransD = cat(1,TransD,sum(A .* ((1:length(A))-1)));
    
    psic_vs_tl = cat(2,psic_vs_tl,psic);
    probl_vs_tl = cat(2,probl_vs_tl,probl);

end
tEnd = toc(tStart);
fprintf('\nElapsed time is %.6f seconds.\n',tEnd)


%%
s.gauge_org_Gs = gauge_org_Gs;
s.matter_org_Gs = matter_org_Gs;
s.psic_vs_tl = psic_vs_tl;
s.probl_vs_tl = probl_vs_tl;
s.density_up_Mt = density_up_Mt;
s.density_dn_Mt = density_dn_Mt;
s.new_obs_LGt = new_obs_LGt;
s.new_obs_LMt = new_obs_LMt;
s.density_Mt_GL = density_Mt_GL;
s.parity_Mt_GL = parity_Mt_GL;
s.G_corr_MtA = G_corr_MtA;
s.gauge_Field_Mt =gauge_Field_Mt;
s.gauge_abs_orderP_Mt = gauge_abs_orderP_Mt;
s.gauge_orderP_Mt = gauge_orderP_Mt;
s.GIPR = GIPR;
s.TransD = TransD;


%% save result
% save([save_path,filesep,...
%     sprintf('analy_deconf_delta@%02dHz.mat',staG0)],...
%     'J','U','tl','staG0','Tilt_E','BDC',...
%     'ham_elems','basis','T','nt',...
%     'gauge_org_Gs',...
%     'matter_org_Gs','psic_vs_tl','probl_vs_tl',...
%     'density_up_Mt','density_dn_Mt','new_obs_LGt',...
%     'new_obs_LMt','density_Mt_GL','parity_Mt_GL',...
%     'G_corr_MtA','gauge_Field_Mt','gauge_abs_orderP_Mt',...
%     'gauge_orderP_Mt','GIPR','TransD')


