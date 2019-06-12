clear; close all;

%% Settings
lam = 1i*sin(linspace(0,1,6400)*pi/2);
lammid = 1i*sin(linspace(0,1,8250)*pi/2);
lamdense = 1i*sin(linspace(0,1,9600)*pi/2);
saveconfigs = 1;
do_generate = 1;
plot_stab_domains = 1;
plot_isb_vs_cores = 1;


%% Generate/Load the optimized stability domains
if do_generate
    configs = cell(1,1); i=0;
    i=i+1; configs{i} = struct('p',4,'N',6, 'nds',[2 6],'lam',lam);
    i=i+1; configs{i} = struct('p',4,'N',10,'nds',[2 10],'lam',lam);
    i=i+1; configs{i} = struct('p',4,'N',14,'nds',[2 14],'lam',lam);
    i=i+1; configs{i} = struct('p',4,'N',18,'nds',[2 18],'lam',lam);
    i=i+1; configs{i} = struct('p',4,'N',22,'nds',[2 8],'lam',lam);
    i=i+1; configs{i} = struct('p',4,'N',26,'nds',[2 26],'lam',lam);
    i=i+1; configs{i} = struct('p',4,'N',28,'nds',[2 28],'lam',lam);

    i=i+1; configs{i} = struct('p',8,'nds',[2 16 18 20],'isbn',.5799);
    i=i+1; configs{i} = struct('p',8,'N',14,'nds',[2 6 8 12],'lam',lam);
    i=i+1; configs{i} = struct('p',8,'N',18,'nds',[2 10 12 14],'lam',lam);
    i=i+1; configs{i} = struct('p',8,'N',22,'nds',[2 4 6 10],'lam',lam);
    i=i+1; configs{i} = struct('p',8,'N',26,'nds',[2 4 8 26],'lam',lam);
    i=i+1; configs{i} = struct('p',8,'N',30,'nds',[2 4 10 30],'lam',lam);
    i=i+1; configs{i} = struct('p',8,'N',34,'nds',[2 4 6 24],'lam',lam);
    i=i+1; configs{i} = struct('p',8,'N',38,'nds',[2 4 10 36],'lam',lam);
    i=i+1; configs{i} = struct('p',8,'N',40,'nds',[2 36:2:40],'lam',lam);

    i=i+1; configs{i} = struct('p',12,'nds',[2 8 12 14 16 20],'isbn',.4515);
%     i=i+1; configs{i} = struct('p',12,'N',14,'nds',[2 4 6 8 10 12]);  % Fully determined system matches here
    i=i+1; configs{i} = struct('p',12,'N',18,'nds',[2 4 6 10 14 16],'lam',lam);
    i=i+1; configs{i} = struct('p',12,'N',22,'nds',[2 4 6 8 12 14],'lam',lam);
    i=i+1; configs{i} = struct('p',12,'N',26,'nds',[2 4 6 8 16 22],'lam',lam);
    i=i+1; configs{i} = struct('p',12,'N',30,'nds',[2 8 10 16 24 26],'lam',lam);
    i=i+1; configs{i} = struct('p',12,'N',34,'nds',[2 8 10 24 26 30],'lam',lam);
    i=i+1; configs{i} = struct('p',12,'N',36,'nds',[2 28:2:36],'lam',lam);

    i=i+1; configs{i} = struct('p',16,'nds',[2 8 10 12 14 16 18 22],'isbn',.4162);
%     i=i+1; configs{i} = struct('p',16,'N',18,'nds',[2 4 6 8 10 12 14 18]);  % Fully determined system wins here
    i=i+1; configs{i} = struct('p',16,'N',22,'nds',[2 4 6 8 10 12 14 18],'lam',lam);
    i=i+1; configs{i} = struct('p',16,'N',26,'nds',[2 4 6 8 10 12 16 20],'lam',lam);
    i=i+1; configs{i} = struct('p',16,'N',30,'nds',[2 18 20 22 24 26 28 30],'lam',lammid);
    i=i+1; configs{i} = struct('p',16,'N',32,'nds',[2 20:2:32],'lam',lamdense);
    
    for i=1:length(configs)
        if isfield(configs{i},'isbn')
            configs{i}.N = max(configs{i}.nds);
            configs{i}.nfs = [];
            configs{i}.optpol = [];            
            configs{i}.ncores = corecount(configs{i}.nds);
        else
            N = configs{i}.N;
            threads = 2:2:N;
            configs{i}.nfs = setdiff(threads,configs{i}.nds);
            configs{i}.optpol = [];
            configs{i}.isbn = 0;
            configs{i}.ncores = corecount(threads);
        end
    end

    parfor i=1:length(configs)
        order=configs{i}.p;
        nds = configs{i}.nds;
        nfs = configs{i}.nfs;    
        maxthread = configs{i}.N;

        fprintf('%d: Order: %d, Maxthread: %d\n', i, order, maxthread);

        if configs{i}.isbn ~= 0
            % Fully determined: compute the coefficints
            y = gbspolys(maxthread);
            [optpol,cd] = rextrap(y(nds,:),nds.^2);
            configs{i}.optpol = fliplr(optpol);
            configs{i}.cd = cd;
            configs{i}.cf = [];
        else
            if order<12
                optfun = @optextrap;
            else
                optfun = @optextrap_vpa;
            end
            try
                [cf,cd,isb,optpol]=optfun(configs{i}.lam,nfs,nds,'tol_bisect',1e-4,'tol_feasible',0);
                isbn = isb/(maxthread+1);
                configs{i}.optpol = optpol;
                configs{i}.isbn = isbn;
                configs{i}.cf = cf;
                configs{i}.cd = cd;
            catch
                warning('error encountered running config %d',i);
            end
        end
    end

    if saveconfigs
        denseconfigs = configs;
        save('denseconfigs.mat','denseconfigs');
    end
else
    load('denseconfigs.mat');
    configs = denseconfigs;
end

%% Plot Stability Domains
if plot_stab_domains
    figyscale = 1.4;
    annotyoffset = -0.018;
    % RK4 Stability Domain
    plot_rk4 = 1;
    if plot_rk4
        npoints = 512;
        [xx,yy] = meshgrid(linspace(-1.05,0.1,npoints),linspace(-1.05,1.05,npoints));
        xi = complex(xx,yy);

        % Display RK4 stability domain
        pol = 1./factorial(0:4);
        u2 = polyval(fliplr(pol),4*xi);
        ua = abs(u2);

        fig = figure(4);
        set(fig,'Position',[100, 100, 800, figyscale*500]);
        ax = getax(fig);
        maketight(ax,.995);
        contour(ax,xx,yy,ua,[1 1],'LineColor',[0.8500, 0.3250, 0.0980],'LineWidth',1.8,'LineStyle','--');
        hold on;
        axlims = [-.8 .1 -1.05 1.05];
        axis(axlims);

        set(ax,'FontSize',24);
        h = xlabel('Re(\xi)'); h.FontSize = 26;
        h = ylabel('Im(\xi)'); h.FontSize = 26;

        xlimits = [-.8 .1];
        ylimits = [-1.05 1.05];
        set(ax,'xtick',xlimits(1):.1:0);
        maketight(ax,.995);

        pos = get(ax,'Position');
        makedim = @(x,y)[(x - min(xlimits))/diff(xlimits) * pos(3) + pos(1), ...
                     (y - min(ylimits))/diff(ylimits) * pos(4) + pos(2), ...
                     0.068 0.065];
        makeannot = @(x,y,str)annotation('textbox',makedim(x,y),'String',str,...
                     'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                     'EdgeColor','none','FaceAlpha',0,'FontSize',24);
        makeannot(-.695,annotyoffset,sprintf('RK'));
    end

    % RK8 Stability Domain
    plot_rk8 = 1;
    if plot_rk8
        npoints = 512;
        [xx,yy] = meshgrid(linspace(-1.05,0.1,npoints),linspace(-1.05,1.05,npoints));
        xi = complex(xx,yy);

        % Display RK8 stability domain
        pol = fliplr(sym2poly(rk8s(1,sym('z'))));
        u2 = polyval(fliplr(pol),13*xi);
        ua = abs(u2);

        fig = figure(8);
        set(fig,'Position',[100, 100, 800, figyscale*500]);
        ax = getax(fig);
        maketight(ax,.995);
        contour(ax,xx,yy,ua,[1 1],'LineColor',[0.8500, 0.3250, 0.0980],'LineWidth',1.8,'LineStyle','--');     
        hold on;
        axlims = [-.5 .05 -1.05 1.05];
        axis(axlims);

        set(ax,'FontSize',24);
        h = xlabel('Re(\xi)'); h.FontSize = 26;
        h = ylabel('Im(\xi)'); h.FontSize = 26;
        
        xlimits = [-.5 .1];
        ylimits = [-1.05 1.05];
        set(ax,'xtick',xlimits(1):.1:0);
        maketight(ax,.995);

        pos = get(ax,'Position');
        makedim = @(x,y)[(x - min(xlimits))/diff(xlimits) * pos(3) + pos(1), ...
                     (y - min(ylimits))/diff(ylimits) * pos(4) + pos(2), ...
                     0.068 0.065];
        makeannot = @(x,y,str)annotation('textbox',makedim(x,y),'String',str,...
                     'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                     'EdgeColor','none','FaceAlpha',0,'FontSize',24);
        makeannot(-.43,annotyoffset,sprintf('RK'));
    end


    % GBS Stability Domains
    haslegend = zeros(4,1);
    for i=1:length(configs)
        order=configs{i}.p;
        maxthread = configs{i}.N;
        optpol = configs{i}.optpol;
        isbn = configs{i}.isbn;
        ncores = configs{i}.ncores;
        if isempty(configs{i}.nfs)
            continue
        end
        if order==4
            indoffset = 0;
            xlimits = [-.8 .1];
            ylimits = [-1.05 1.05];
            xints = [-.761 -.580 -.479  -.414  -.370  -.328 -.286];
            xoffs = [-.01/7 0.0    .005/7 .01/3  .02/3 .03/3  .04/2];
            xints = xints+xoffs;
        elseif order==8
            indoffset = 8;
            xlimits = [-.5 .1];
            ylimits = [-1.05 1.05];
            xints = [-.458 -.375 -.353 -.317 -.294 -.273 -.256 -.226];
        elseif order==12
            indoffset = 16;
            xlimits = [-.5 .05];
            ylimits = [-1.05 1.05];
            xints = [-.484 -.4432 -.403 -.367 -.343 -.319 -.291];
            xints = xints+.0035;
        elseif order==16
            indoffset = 23;
            xlimits = [-.5 .05];
            ylimits = [-1.05 1.05];
            xints = [-.457 -.437 -.402 -.383 -.348];
            xints = xints+.0035;
        end
        xints = xints-.002;
        axlims = [xlimits ylimits];
        ndsstr = regexprep(num2str(configs{i}.nds),' +','_');
        legendtext = sprintf('  %dth Order',order);
        filename = sprintf('plots/reduced/stab_domain__p_%d__N_%d__nds_%s.png',order,maxthread,ndsstr);

        fig = figure(order);
        set(fig,'Position',[100, 100, 800, figyscale*500]);
        ax = getax(fig);
        maketight(ax,.995);
        plot_stability_domain(ax,optpol,[],0);
        axis(axlims);
        set(ax,'FontSize',24);
        set(ax,'xtick',xlimits(1):.1:0);
        h = xlabel('Re(\xi)'); h.FontSize = 26;
        h = ylabel('Im(\xi)'); h.FontSize = 26;
        maketight(ax,.995);
        if ~haslegend(order/4)
            pos = get(ax,'Position');
            lgd = annotation('textbox',[pos(1)+.025*pos(3),pos(2)+.89*pos(4),.068,.065],'String',legendtext,...
                             'FitBoxToText','on','BackgroundColor',[1 1 1],'FaceAlpha',1);
            lgd.FontSize = 26;
            haslegend(order/4) = 1;
        end


        axis([xlimits ylimits]); 
        pos = get(ax,'Position');
        makedim = @(x,y)[(x - min(xlimits))/diff(xlimits) * pos(3) + pos(1), ...
                         (y - min(ylimits))/diff(ylimits) * pos(4) + pos(2), ...
                         0.068 0.065];
        makeannot = @(x,y,str)annotation('textbox',makedim(x,y),'String',str,...
                         'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                         'EdgeColor','none','FaceAlpha',0,'FontSize',24);
        makeannot(xints(i-indoffset),annotyoffset,sprintf('%d',ncores));
    end

    orders = [4 8 12 16];
    for i=1:length(orders)
        order = orders(i); 
        fig = figure(order);
        ax = getax(fig);
        l1 = plot(ax,[0 0],[-1.05 1.05],'k--','LineWidth',1.0,'HandleVisibility','off');
        uistack(l1,'bottom');
        l2 = plot(ax,[-1. .1],[0 0],'k--','LineWidth',0.8,'HandleVisibility','off');
        uistack(l2,'bottom');

        maketight(ax,.995);
        set(fig,'Color','w');

        paper_export_fig(ax,sprintf('stabdomains_p_%d',order));
    end
end

%% ISB vs Core Curves
if plot_isb_vs_cores
    fig = figure;
    ax = getax(fig);
    xlimits = [-.2 12.4];
    ylimits = [.2 1];
    pos = get(ax,'Position');
    makedim = @(x,y)[(x - min(xlimits))/diff(xlimits) * pos(3) + pos(1) + .004, ...
                     (y - min(ylimits))/diff(ylimits) * pos(4) + pos(2) - .035, ...
                     0.0678571428571435 0.0642857142857143];
    makeannot = @(x,y,str)annotation('textbox',makedim(x,y),'String',str,...
                     'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                     'EdgeColor','none','FaceAlpha',0,'FontSize',12);
    orders = [4 8 12 16];
    markers = {':v',':s',':^',':o'};
    colors = linecolors();
    linewidth = 1.4;
    for j=1:length(orders)
        order = orders(j);
        cores = [];
        isbs = [];
        for i=1:length(configs)
            if configs{i}.p==order
                cores = [cores configs{i}.ncores];
                isbs = [isbs configs{i}.isbn];
            end
        end
        plot(cores,isbs,markers{j},'LineWidth',linewidth,'Color',colors(j,:)); hold on;
        axis([xlimits ylimits]);
        makeannot(cores(end),isbs(end),sprintf('GBS_{%d}',order));
    end
    plot(1,sqrt(2)/2,'*','LineWidth',linewidth,'Color',colors(j+1,:));  % RK4
    plot(1,.2848,'x','LineWidth',linewidth,'Color',colors(j+2,:));      % RK8
    axis([xlimits ylimits]);
    makeannot(1,sqrt(2)/2,'RK_4');
    makeannot(1,.2848,'RK_8');

    grid on;
    ax.FontSize = 10;
    h = xlabel('Number of Cores'); h.FontSize = 12;
    h = ylabel('ISB_n'); h.FontSize = 12;
    set(ax,'xtick',1:11)
    set(fig,'Color','w');

    paper_export_fig(ax,'isb_vs_cores',1);
end

function ax = getax(fig)
    ax = findall(fig,'type','axes');
    if length(ax)>1
        ax = ax(1);
    end
    if isempty(ax)
        ax = axes(fig);
    end
end

