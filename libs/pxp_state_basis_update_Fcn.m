function [basis] = pxp_state_basis_update_Fcn(st_up,st_dn,basis,st_orig)
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
if nargin < 4
    st_orig = ones(size(st_dn));
end


%% search the index of the initial state.
nMax = basis.nMax;
idxlt = basis.idxstatel;

%
idxIntl_up = base2dec(strrep(num2str(st_up),' ',''),nMax+1);
idxIntl_dn = base2dec(strrep(num2str(st_dn),' ',''),nMax+1);

[~,idx_st] = ismember([idxIntl_up,idxIntl_dn],idxlt,'rows');


%%
stateA = basis.state_dn;
st_pxp_init = st_dn;

%
n_bs = basis.n_bs;
L = basis.L;


%% initial gauge field, matter field, and Gauss law
% initial gauge field, in given lattice sites
gauge_init = zeros(1,L);

gauge_init(1:2:end) = 0.5 * ((st_pxp_init(1:2:end) == 2) - 0.5)*2;
gauge_init(2:2:end) = -0.5 * ((st_pxp_init(2:2:end) == 2) - 0.5)*2;

% initial matter field, the start and end matter sites are defined as zero
% at beginning
matter_init = zeros(1,L+1);
matter_init(2:end-1) = gauge_init(2:end) - gauge_init(1:end-1);

% the edge of the guage field, these will keep as fixed value 
gauge_edge = zeros(1,2);
gauge_edge(1) = gauge_init(1);
gauge_edge(2) = gauge_init(end);


%% convert original basis into gauge field under PXP model
gauge_MtA = nan(size(stateA));

stateA_OD = stateA(:,1:2:end);
stateA_EN = stateA(:,2:2:end);

gauge_MtA_OD = nan(size(stateA_OD));
gauge_MtA_EN = nan(size(stateA_EN));

gauge_MtA_OD(stateA_OD == 2) = 0.5;
gauge_MtA_OD((stateA_OD == 0) | (stateA_OD == 1)) = -0.5;
gauge_MtA_OD(isnan(sum(gauge_MtA_OD,2)),:) = NaN;

gauge_MtA_EN(stateA_EN == 2) = -0.5;
gauge_MtA_EN((stateA_EN == 0) | (stateA_EN == 1)) = 0.5;
gauge_MtA_EN(isnan(sum(gauge_MtA_EN,2)),:) = NaN;

gauge_MtA(:,1:2:end) = gauge_MtA_OD;
gauge_MtA(:,2:2:end) = gauge_MtA_EN;
gauge_MtA(isnan(sum(gauge_MtA,2)),:) = NaN;

% idx_gauge_nonNaN = ~isnan(sum(gauge_MtA,2));
% sum(idx_gauge_nonNaN)

% then the full guage field
gauge_All = [repmat(gauge_edge(1),n_bs,1),...
    gauge_MtA,repmat(gauge_edge(2),n_bs,1)];


%% for matter field
matter_All = gauge_All(:,2:end) - gauge_All(:,1:end-1);
matter_MtA = matter_All(:,2:end-1);
% idx_matter_nonNaN = ~isnan(sum(matter_MtA,2));
% sum(idx_matter_nonNaN)


%% gauss law for initial state
gl_init = st_pxp_init * (1:length(st_pxp_init))' ...
    + 0.5 * st_pxp_init * (st_pxp_init'-1);
gl_stateA = stateA * (1:length(st_pxp_init))' ...
    + 0.5 * sum(stateA.*(stateA-1),2);
gl_idx_init_lt = gl_stateA == gl_init;


%% index search part
idx_gauge = ~isnan(sum(gauge_MtA,2));

matter_od_All = matter_All(:,1:2:end);
matter_en_All = matter_All(:,2:2:end);

nm_od = size(matter_od_All,2);
nm_en = size(matter_en_All,2);

idx_od = sum((matter_od_All==0) + (matter_od_All==-1),2) == nm_od;
idx_en = sum((matter_en_All==0) + (matter_en_All==1),2) == nm_en;

idx_all = idx_gauge & idx_od & idx_en;
GL_idx_lt = gl_idx_init_lt & idx_all;


%% original state
st_pxp_orig = st_orig;
gauge_st_pxp_orig = zeros(size(st_pxp_orig));

st_pxp_orig_OD = st_orig(1:2:end);
st_pxp_orig_EN = st_orig(2:2:end);

gauge_st_pxp_orig_OD = zeros(size(st_pxp_orig_OD));
gauge_st_pxp_orig_EN = zeros(size(st_pxp_orig_EN));

gauge_st_pxp_orig_OD(st_pxp_orig_OD == 2) = 0.5;
gauge_st_pxp_orig_OD((st_pxp_orig_OD == 0) | (st_pxp_orig_OD == 1)) = -0.5;

gauge_st_pxp_orig_EN(st_pxp_orig_EN == 2) = -0.5;
gauge_st_pxp_orig_EN((st_pxp_orig_EN == 0) | (st_pxp_orig_EN == 1)) = 0.5;

gauge_st_pxp_orig(1:2:end) = gauge_st_pxp_orig_OD;
gauge_st_pxp_orig(2:2:end) = gauge_st_pxp_orig_EN;

matter_st_pxp_orig = gauge_st_pxp_orig(2:end) - gauge_st_pxp_orig(1:end-1);

gauge_st_pxp_origA = [gauge_edge(1),gauge_st_pxp_orig,gauge_edge(2)];
matter_st_pxp_origA = [0,matter_st_pxp_orig,0];


%% return value
basis.idx_init_st = idx_st;

basis.gauge_All = gauge_All;
basis.matter_All = matter_All;

basis.matter_init = matter_init;
basis.gauge_edge = gauge_edge;

basis.gauge_MtA = gauge_MtA;
basis.matter_MtA = matter_MtA;

basis.GL_idx_lt = GL_idx_lt;

basis.gauge_st_pxp_orig = gauge_st_pxp_orig;
basis.matter_st_pxp_orig = matter_st_pxp_orig;
basis.gauge_st_pxp_origA = gauge_st_pxp_origA;
basis.matter_st_pxp_origA = matter_st_pxp_origA;

basis.gauge_st_pxp_init = gauge_st_pxp_orig;
basis.matter_st_pxp_init = matter_st_pxp_orig;

