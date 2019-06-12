clear; close all;

orders = [4 8];
isbns = [1/sqrt(2) 0.2848];
pol_rk4 = 1./factorial(0:4);
pol_rk8 = fliplr(sym2poly(rk8s(1,sym('z'))));

y = zeros(2,length(pol_rk8));
y(1,1:length(pol_rk4)) = pol_rk4;
y(2,:) = pol_rk8;

figure;
fig = gcf;
set(fig,'Color',[1,1,1]);
set(fig,'Position',[100 100 1000 600]);

nplots = size(y,1);
for i=1:nplots
    xlimits = [-.8 .1];
    ylimits = [-1.05 1.05];

    ax = subplot(1,2,i);
    plot_stability_domain(ax,fliplr(y(i,:)));

    axis(ax,[xlimits ylimits]);
    ax.FontSize = 16;
    h=xlabel('Re(\xi)'); h.FontSize = 18;
    h=ylabel('Im(\xi)'); h.FontSize = 18;

    order = orders(i);
    isbn = isbns(i);
    isbprint = floor(isbn*1e4)/1e4;
    legendtext = sprintf('  RK_{%d},   ISB_n = %1.4f',order,isbprint);

    pos = get(ax,'Position');
    annotation('textbox',[pos(1)+.025*pos(3),pos(2)+.89*pos(4),.068,.065],'String',legendtext,...
                     'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                     'EdgeColor',[0 0 0],'FaceAlpha',1,'FontSize',18);
end
paper_export_fig(fig,'stabdomain_rk_4_8');

