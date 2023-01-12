function [s_corr] = correlation_statistics_fcn(psic,basis)

%%
s_corr = struct();


%%
% n_bs = basis.n_bs;
L = basis.L;
nMax = basis.nMax;
idxstatel = basis.idxstatel;

state_up = basis.state_up;
state_dn = basis.state_dn;


%%
% selection rules current
idx_sum_eq2 = (state_up + state_dn) == 2;

% correlation, 
% S_i^+ * S_j^- = A_i^+ * B_i * A_j * B_j^+
% s_corr = struct();


%%

for ss = 0:L-1
    num_corr_cur = L-ss;
    corr_lt_cur_raw = zeros(1,num_corr_cur);
    corr_lt_cur_sel = zeros(1,num_corr_cur);
    corr_cur_strings = cell(1,num_corr_cur);
    field_cur = ['corr_dis_',num2str(ss)];
    
    for kk = 1:L-ss
        ia = kk;
        ib = ia + ss;
        field_sub = ['SpSm_',num2str(ia),'_',num2str(ib)];
        corr_cur_strings{kk} = field_sub;
        
        t_jCur = tic;
        fprintf('\nstatistics S+S- between site %02d and %02d.\n',ia,ib)

        if ss == 0
            stat_sub_up_raw = state_up(:,ia);
            stat_sub_dn_raw = state_dn(:,ia);
            
            % A_ia^+ * B_ia * A_ib * B_ib^+ = n_A * (n_B+1)
            idx_up_raw = stat_sub_up_raw > 0;
            idx_dn_raw = stat_sub_dn_raw < nMax;
            idx_raw = idx_up_raw & idx_dn_raw;
            
            corr_coeff_cur_all = stat_sub_up_raw .* (stat_sub_dn_raw+1);
            
            % raw case
            psi_bra_raw = conj(psic(idx_raw));
            psi_ket_raw = psic(idx_raw);
            corr_coeff_raw = corr_coeff_cur_all(idx_raw);
            corr_raw = psi_bra_raw.' * (corr_coeff_raw.*psi_ket_raw);
            
            % selected case
            idx_keep = idx_sum_eq2(:,ia);
            psic_sel = psic .* idx_keep;
            psic_sel = psic_sel/sqrt(psic_sel'*psic_sel);
            idx_sel = idx_raw & idx_keep;
            psi_bra_sel = conj(psic_sel(idx_sel));
            psi_ket_sel = psic_sel(idx_sel);
            corr_coeff_sel = corr_coeff_cur_all(idx_sel);
            corr_sel = psi_bra_sel.' * (corr_coeff_sel.*psi_ket_sel);
            
            corr_lt_cur_raw(kk) = corr_raw;
            corr_lt_cur_sel(kk) = corr_sel;
            
        else
            stat_sub_up_raw = state_up(:,[ia,ib]);
            stat_sub_dn_raw = state_dn(:,[ia,ib]);

            % A_ia^+ * B_ia * A_ib * B_ib^+
            idx_up_ia_raw = stat_sub_up_raw(:,1) < nMax;
            idx_dn_ia_raw = stat_sub_dn_raw(:,1) > 0;
            idx_up_ib_raw = stat_sub_up_raw(:,2) > 0;
            idx_dn_ib_raw = stat_sub_dn_raw(:,2) < nMax;
            idx_raw = idx_up_ia_raw & idx_dn_ia_raw & ...
                idx_up_ib_raw & idx_dn_ib_raw;

            % 01: raw case
            psi_ket_raw = psic(idx_raw);
            
            stat_sub_up = stat_sub_up_raw(idx_raw,:);
            stat_sub_dn = stat_sub_dn_raw(idx_raw,:);
            corr_coeff_raw = sqrt(stat_sub_up(:,1)+1) ...
                .* sqrt(stat_sub_dn(:,1)) .* sqrt(stat_sub_up(:,2)) ...
                .* sqrt(stat_sub_dn(:,2)+1);

            %
            state_up_nxt_raw = state_up(idx_raw,:);
            state_dn_nxt_raw = state_dn(idx_raw,:);
            state_up_nxt_raw(:,ia) = state_up_nxt_raw(:,ia) + 1;
            state_up_nxt_raw(:,ib) = state_up_nxt_raw(:,ib) - 1;
            state_dn_nxt_raw(:,ia) = state_dn_nxt_raw(:,ia) - 1;
            state_dn_nxt_raw(:,ib) = state_dn_nxt_raw(:,ib) + 1;
            idx_lt_up_nxt_raw = state_up_nxt_raw * ((nMax+1).^(L-1:-1:0))'; 
            idx_lt_dn_nxt_raw = state_dn_nxt_raw * ((nMax+1).^(L-1:-1:0))';
            [~,idx_nxt_raw] = ismember([idx_lt_up_nxt_raw,idx_lt_dn_nxt_raw],...
                idxstatel,'rows');
            psi_bra_raw = psic(idx_nxt_raw);
            
            corr_raw = psi_bra_raw' * (corr_coeff_raw.*psi_ket_raw);
            
            % 02: select case
            idx_keep = idx_sum_eq2(:,ia) & idx_sum_eq2(:,ib);
            psic_sel = psic .* idx_keep;
            psic_sel = psic_sel/sqrt(psic_sel'*psic_sel);
            idx_sel = idx_raw & idx_keep;
            psi_in_sel = psic_sel(idx_sel);
            
            stat_sub_up = stat_sub_up_raw(idx_sel,:);
            stat_sub_dn = stat_sub_dn_raw(idx_sel,:);
            corr_coeff_sel = sqrt(stat_sub_up(:,1)+1) ...
                .* sqrt(stat_sub_dn(:,1)) .* sqrt(stat_sub_up(:,2)) ...
                .* sqrt(stat_sub_dn(:,2)+1);

            %
            state_up_nxt_sel = state_up(idx_sel,:);
            state_dn_nxt_sel = state_dn(idx_sel,:);
            state_up_nxt_sel(:,ia) = state_up_nxt_sel(:,ia) + 1;
            state_up_nxt_sel(:,ib) = state_up_nxt_sel(:,ib) - 1;
            state_dn_nxt_sel(:,ia) = state_dn_nxt_sel(:,ia) - 1;
            state_dn_nxt_sel(:,ib) = state_dn_nxt_sel(:,ib) + 1;
            idx_lt_up_nxt_sel = state_up_nxt_sel * ((nMax+1).^(L-1:-1:0))'; 
            idx_lt_dn_nxt_sel = state_dn_nxt_sel * ((nMax+1).^(L-1:-1:0))';
            [~,idx_nxt_sel] = ismember([idx_lt_up_nxt_sel,idx_lt_dn_nxt_sel],...
                idxstatel,'rows');
            psi_out_sel = conj(psic_sel(idx_nxt_sel));
            
            corr_sel = psi_out_sel.' * (corr_coeff_sel.*psi_in_sel);
            
            corr_lt_cur_raw(kk) = corr_raw;
            corr_lt_cur_sel(kk) = corr_sel;
            
        end
        
        s_corr.(field_cur).(field_sub).corr_raw = corr_raw;
        s_corr.(field_cur).(field_sub).corr_sel = corr_sel;
        
        tDJc = toc(t_jCur);
        fprintf('elapsed time is %.6f seconds.\n',tDJc)

    end

    s_corr.(field_cur).num_corr_cur = num_corr_cur;
    s_corr.(field_cur).corr_lt_cur_raw = corr_lt_cur_raw;
    s_corr.(field_cur).corr_lt_cur_sel = corr_lt_cur_sel;
    s_corr.(field_cur).corr_cur_strings = corr_cur_strings;
    
end


% s_corr.corr_dis_0.SpSm_1_1
