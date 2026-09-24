% 示例：生成一个含NaN的随机数据
% A = rand(30,40);
% A(10:20, 15:25) = NaN;          % 制造一个空洞区域
% A(5:25, 5:10) = NaN;            % 制造另一个区域
j = 10;
filepat = 'accu_planar_0709_xmoon\';
filepat1 = strcat(filepat,'J',num2str(J0(j)),'_r2',num2str(r2),'_z0',num2str(z0),...
        '_dz0',num2str(dz0));
matstr = strcat(filepat1,'\','LD_',num2str(nonlinear_id),'_',LD_str,...
'_expalp',num2str(exp_alp),'_changeT_','0709','.mat');
load(matstr,'LD_FTLE_all','LD_a2_all','Np_prl_all_t',...
        'KEarth_all_t','KMoon_all_t','KEscape_all_t');
% load()
%%
kt = 5;
LD_a2_mesh=LD_a2_all(:,:,kt) ;
Np_prl_all=Np_prl_all_t(:,:,kt) ;
KEarth_all=KEarth_all_t(:,:,kt);
KMoon_all=KMoon_all_t(:,:,kt);
KEscape_all = KEscape_all_t(:,:,kt);
KXMoon_all = KXMoon_all_t(:,:,kt);

LD_FTLE_mesh=LD_FTLE_all(:,:,kt);

LD_FTLE_mesh2 = LD_FTLE_mesh;
LD_FTLE_mesh2(LD_a2_mesh<=0.55) = nan;
% 提取边界
BW = ~isnan(LD_FTLE_mesh2);                 % 二值图像
boundaries = bwboundaries(BW, 'noholes');  % 'noholes' 只返回外边界，忽略内部空洞

% 绘图查看
figure;
% imshow(BW, 'InitialMagnification', 'fit'); hold on;
% for k = 1:length(boundaries)
%     b = boundaries{k};
%     plot(b(:,2), b(:,1), 'r', 'LineWidth', 2);  % 注意坐标：列(x) 行(y)
% end
% title('非NaN区域的边界');

% fig_LD = figure('Position', [46.60,192.2,949.60,560], 'Color', 'white');
hAx = axes;
LD_field=surf(hAx,alp_mesh/rad,Psi_mesh/rad,LD_FTLE_mesh2, 'EdgeColor', 'none');
colormap(hAx,jet);colorbar(hAx);hold on;
xlim([min(alp/rad),max(alp/rad)]); ylim([min(Psi/rad),max(Psi/rad)]);
for k = 1:length(boundaries)
    b = boundaries{k};
    plot3(alp(b(:,2))/rad,Psi(b(:,1))/rad,8*ones(size(b(:,1))), 'r', 'LineWidth', 2.8);  % 注意坐标：列(x) 行(y)
end
xlabel('\alpha \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
ylabel('\Psi \rm[deg]','fontsize',25,'fontname','Times New Roman','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');

