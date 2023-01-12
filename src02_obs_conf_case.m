close all
clc

fileAll = dir([pwd,filesep,'analy_res_conf',filesep,'analy_*.mat']);
nfls = numel(fileAll);


%%
delta_lt = [];

new_obs_Gt = [];
new_obs_Mt = [];

G_Mt = [];
T_Mt = [];

E_Mt = [];
tl = [];

for kk = 1:nfls
    fname = fileAll(kk).name;
    dname = fileAll(kk).folder;
    
    delta_lt = cat(1,delta_lt,...
        str2double(fname(end-7:end-6)));
    
    t_matC = load([dname,filesep,fname]);
    new_obs_Gt = cat(2,new_obs_Gt,t_matC.new_obs_LGt);
    new_obs_Mt = cat(2,new_obs_Mt,t_matC.new_obs_LMt);
    
    tl = t_matC.tl;
    
end


%%
nC = length(delta_lt);


% Color define
A = {'264653','2a9d8f','e9c46a','f4a261','e76f51'};

colorSE = zeros(numel(A),3);
for kk = 1:numel(A)
    colorSE(kk,:) = hex2rgb(A{kk});
end

[XE,YE] = meshgrid(1:3,1:numel(A));
[XqE,YqE] = meshgrid(1:3,linspace(1,numel(A),nC));
% Interpolation for 1-D gridded data
VqE = interp2(XE,YE,colorSE,XqE,YqE,'spline');
VqE(VqE<0) = 0;
VqE(VqE>1) = 1;



%%
figure('Color','w','Position',[120 120 422 320])
hold on
for kk = 1:nfls
    plot(tl*1000,new_obs_Gt(:,kk),'-',...
        'Color',[0 0.4470 0.7410],'LineWidth',3,...
        'DisplayName',sprintf('\\delta = %02d Hz',delta_lt(kk)))
end
hold off
box on

lgd = legend('Show');
lgd.Location = 'best';

lgd.BoxFace.ColorType='truecoloralpha';
lgd.BoxFace.ColorData=uint8(255*[1 1 1 0.10]');

ax = gca;
ax.FontSize = 18;
xlim([-0.5,50])
ylim([0.5,2])
xlabel('Evolution time \itt\rm (ms)','FontSize',20)
ylabel('\langle\itO(t)\rangle','FontSize',20)


export_fig('Numerical_Obs_Ot.pdf')

% export_fig('-r300','Delta_EField_Numeric.png')


%%
figure('Color','w','Position',[420 120 422 320])
hold on
for kk = 1:nfls
    if kk == 1
        plot(tl(1:end)*1000,new_obs_Mt(1:end,kk),'-',...
            'Color',[0 0.4470 0.7410],'LineWidth',2,...
            'DisplayName',sprintf('\\delta = %02d Hz',delta_lt(kk)))
    else
        plot(tl*1000,new_obs_Mt(:,kk),'-',...
            'Color',VqE(kk,:),'LineWidth',2,...
            'DisplayName',sprintf('\\delta = %02d Hz',delta_lt(kk)))
    end
end
hold off
box on

lgd = legend('Show');
lgd.Location = 'best';

lgd.BoxFace.ColorType='truecoloralpha';
lgd.BoxFace.ColorData=uint8(255*[1 1 1 0.10]');

ax = gca;
ax.FontSize = 18;
xlim([-0.5,50])
% ylim([-0.01,0.25])
ylim([1.5,3])
xlabel('Evolution time \itt\rm (ms)','FontSize',20)
ylabel('\langle\itM(t)\rangle','FontSize',20)

export_fig('Numerical_Obs_Mt.pdf')

% export_fig('-r300','Delta_mField_Numeric.png')



%%

ratio_mt = 2*new_obs_Gt./new_obs_Mt;

figure('Color','w','Position',[720 120 422 320])
hold on
for kk = 1:nfls

    plot(tl*1000,ratio_mt(:,kk),'-',...
        'Color',[0 0.4470 0.7410],'LineWidth',2,...
        'DisplayName',sprintf('\\delta = %02d Hz',delta_lt(kk)))

end
hold off
box on

lgd = legend('Show');
lgd.Location = 'best';

lgd.BoxFace.ColorType='truecoloralpha';
lgd.BoxFace.ColorData=uint8(255*[1 1 1 0.10]');

ax = gca;
ax.FontSize = 18;
xlim([-0.5,50])
% ylim([-0.01,0.25])
ylim([0.5,2])
xlabel('Evolution time \itt\rm (ms)','FontSize',20)
ylabel('\langle\it2 O(t) / M(t)\rangle','FontSize',20)

export_fig('Numerical_Obs_ratio.pdf')

% export_fig('-r300','Delta_Ratio_Numeric.png')

return
