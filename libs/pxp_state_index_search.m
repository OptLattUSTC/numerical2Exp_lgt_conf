function [basis] = pxp_state_index_search(st_up,st_dn,basis,st_orig)


%%
if nargin < 4
    st_orig = ones(size(st_dn));
end


%%
nMax = basis.nMax;
idxlt = basis.idxstatel;

st_cur = st_dn;

%%
% idxlt = basis.idxIntl;
stateA = basis.state_dn;



%%
% defination in PXP model
% site occupation vs gauge field
%
%          |   ODD  |   EVEN
%   ------ | ------ | -------
%     ->   |    2   |   0/1
%     <-   |   0/1  |    2
%
% with -> = 0.5; <- = -0.5
%


%
gauge_MtA = nan(size(stateA));
matter_MtA = gauge_MtA(:,1:end-1);

stateA_OD = stateA(:,1:2:end);
st_cur_EN = stateA(:,2:2:end);

gauge_MtA_OD = nan(size(stateA_OD));
gauge_MtA_EN = nan(size(st_cur_EN));
% gauge_MtA_EN = gauge_MtA_OD;

gauge_MtA_OD(stateA_OD == 2) = 0.5;
gauge_MtA_OD((stateA_OD == 0) | (stateA_OD == 1)) = -0.5;
gauge_MtA_OD(isnan(sum(gauge_MtA_OD,2)),:) = NaN;

gauge_MtA_EN(st_cur_EN == 2) = -0.5;
gauge_MtA_EN((st_cur_EN == 0) | (st_cur_EN == 1)) = 0.5;
gauge_MtA_EN(isnan(sum(gauge_MtA_EN,2)),:) = NaN;

gauge_MtA(:,1:2:end) = gauge_MtA_OD;
gauge_MtA(:,2:2:end) = gauge_MtA_EN;
gauge_MtA(isnan(sum(gauge_MtA,2)),:) = NaN;

matter_MtA = gauge_MtA(:,2:end) - gauge_MtA(:,1:end-1);


%%
idxIntl_up = base2dec(strrep(num2str(st_up),' ',''),nMax+1);
idxIntl_dn = base2dec(strrep(num2str(st_dn),' ',''),nMax+1);

[~,idx_st] = ismember([idxIntl_up,idxIntl_dn],idxlt,'rows');


%% gauss law for initial state
GL_init = st_cur * (1:length(st_cur))' ...
    + 0.5 * st_cur * (st_cur'-1);
GL_stateA = stateA * (1:length(st_cur))' ...
    + 0.5 * sum(stateA.*(stateA-1),2);
GL_idx_lt = GL_stateA == GL_init;

%
idx_non_nan = ~isnan(sum(gauge_MtA,2));
GL_idx_lt = GL_idx_lt & idx_non_nan;


%%
% st_cur
gauge_st_cur = zeros(size(st_cur));

st_cur_OD = st_orig(1:2:end);
st_cur_EN = st_orig(2:2:end);

gauge_st_cur_OD = zeros(size(st_cur_OD));
gauge_st_cur_EN = zeros(size(st_cur_EN));
% gauge_st_cur_EN = gauge_st_cur_OD;

gauge_st_cur_OD(st_cur_OD == 2) = 0.5;
gauge_st_cur_OD((st_cur_OD == 0) | (st_cur_OD == 1)) = -0.5;

gauge_st_cur_EN(st_cur_EN == 2) = -0.5;
gauge_st_cur_EN((st_cur_EN == 0) | (st_cur_EN == 1)) = 0.5;

gauge_st_cur(1:2:end) = gauge_st_cur_OD;
gauge_st_cur(2:2:end) = gauge_st_cur_EN;

matter_st_cur = gauge_st_cur(2:end) - gauge_st_cur(1:end-1);


%%
basis.GL_idx_lt = GL_idx_lt;
basis.idx_init_st = idx_st;

basis.gauge_MtA = gauge_MtA;
basis.matter_MtA = matter_MtA;

basis.gauge_st_cur = gauge_st_cur;
basis.matter_st_cur = matter_st_cur;

