function stat_n = pxpM_site_fill_statFcn(psic,basis)


%%
stat_n = struct();


%%
n_bs = basis.n_bs;
L = basis.L;

state_up = basis.state_up;
state_dn = basis.state_dn;


%% probability list
prob_ltC = abs(conj(psic).*psic);


%%
density_up_ltC = sum(state_up .* repmat(prob_ltC,1,L),1);
parity_up_ltC = sum(mod(state_up,2) .* repmat(prob_ltC,1,L),1);

density_dn_ltC = sum(state_dn .* repmat(prob_ltC,1,L),1);
parity_dn_ltC = sum(mod(state_dn,2) .* repmat(prob_ltC,1,L),1);


%%
stat_n.L = L;
stat_n.n_bs = n_bs;

stat_n.prob_ltC = prob_ltC;

stat_n.density_up_ltC = density_up_ltC;
stat_n.parity_up_ltC = parity_up_ltC;

stat_n.density_dn_ltC = density_dn_ltC;
stat_n.parity_dn_ltC = parity_dn_ltC;



%%
% defination in PXP model
% site occupation vs gauge field
%
%          |   ODD  |   EVEN
%   ------ | ------ | -------
%     ->   |    2   |   0/1
%     <-   |   0/1  |    2
%


%%
ns = n_bs;
L = basis.L;
state_mtA = basis.state_dn;
GL_idx_lt = basis.GL_idx_lt;

matter_MtA = basis.matter_MtA;
gauge_MtA = basis.gauge_MtA;


%% probability list
prob_ltC_GL = prob_ltC .* GL_idx_lt;
prob_ltC_GL = prob_ltC_GL./sum(prob_ltC_GL);

density_ltC_GL = sum(state_mtA .* repmat(prob_ltC_GL,1,L),1);
parity_ltC_GL = sum(mod(state_mtA,2) .* repmat(prob_ltC_GL,1,L),1);


%%
stat_n.GL_idx_lt = GL_idx_lt;
stat_n.prob_ltC_GL = prob_ltC_GL;
stat_n.density_ltC_GL = density_ltC_GL;
stat_n.parity_ltC_GL = parity_ltC_GL;


%%
guage_Mt_Sub = gauge_MtA(basis.GL_idx_lt,:);
prob_ltC_Sub = prob_ltC_GL(basis.GL_idx_lt);

gauge_Sub_lt = prob_ltC_Sub' * guage_Mt_Sub;


%%
matter_DensitylC = sum(matter_MtA .* repmat(prob_ltC_GL,1,L-1),1,'omitnan');
gauge_Density_lt = sum(gauge_MtA .* repmat(prob_ltC_GL,1,L),1,'omitnan');

% gauge_DensitylC = sum(abs(gauge_MtA) .* repmat(prob_ltC_GL,1,L),1,'omitnan');
% gauge_DensitylC = sum(abs(gauge_MtA) .* repmat(prob_ltC_GL,1,L),1,'omitnan');

gauge_DensitylC = sum(abs(sum(gauge_MtA,2) .* prob_ltC_GL),'omitnan')/L;

matter_Density_Init = abs(basis.matter_st_cur);
matter_MtA0 = abs(matter_MtA);

% correlation function
G_dis_r = zeros(1,L-1);
for kk = 1:L-1
    r = kk - 1;
    G_cur = zeros(ns,1);
    for jj = 1:L-1
        try
            G_cur = G_cur + (matter_MtA0(:,jj) - matter_Density_Init(jj)) ...
                .* (matter_MtA0(:,jj+r) - matter_Density_Init(jj+r));
        catch
            continue
        end
        
    end
    
    G_dis_r(kk) = sum(G_cur .* prob_ltC_GL,1,'omitnan')/(L-kk);
    
end

stat_n.matter_DensitylC = matter_DensitylC;
stat_n.G_dis_r = G_dis_r;
stat_n.gauge_DensitylC = gauge_DensitylC;

stat_n.gauge_Density_lt = gauge_Density_lt;
stat_n.gauge_Sub_lt = gauge_Sub_lt;


