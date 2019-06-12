clear; close all;

%% Script parameters
saveconfigs = 1;
do_generate = 1;
do_rationalize = 1;
ndigits = 256;
digits(ndigits);

%% Construct the optimization curve
npoints = 6400;

%% Construct the configurations
configs = cell(1,1); i=0;
i=i+1; configs{i} = struct('p',8, 'N',22,'nds',[2 4 6 10],'numapprox',3, ...
                           'tol',0,'lam',make_curve(npoints,.5,.1,.01,1e-3));      % 6 cores
i=i+1; configs{i} = struct('p',8, 'N',30,'nds',[2 26:2:30],'numapprox',2, ...
                           'tol',1e-12,'lam',make_curve(npoints,.4,.2,.01,1e-3));  % 8 cores
i=i+1; configs{i} = struct('p',12,'N',30,'nds',[2 8 10 16 24 26],'numapprox',3, ...
                           'tol',0,'lam',make_curve(npoints,.5,.2,.01,4e-4));      % 8 cores

for i=1:length(configs)
    N = configs{i}.N;
    threads = 2:2:N;
    configs{i}.nfs = setdiff(threads,configs{i}.nds);
    configs{i}.optpol = [];
    configs{i}.isbn = 0;
    configs{i}.ncores = corecount(threads);
end

%% Run the optimizer
if do_generate
    parfor i=1:length(configs)
        digits(ndigits);
        order=configs{i}.p;
        nds = configs{i}.nds;
        nfs = configs{i}.nfs;    
        maxthread = configs{i}.N;
        lam = configs{i}.lam;
        tol_feasible = configs{i}.tol;

        fprintf('%d: Order: %d, Maxthread: %d\n', i, order, maxthread);

        if order<12
            optfun = @optextrap;
        else
            optfun = @optextrap_vpa;
        end
        try
            [cf,cd,isb,optpol]=optfun(lam,nfs,nds,'tol_bisect',1e-3,'tol_feasible',tol_feasible);
            isbn = isb/(maxthread+1);
            configs{i}.optpol = optpol;
            configs{i}.isbn = isbn;
            configs{i}.cf = cf;
            configs{i}.cd = cd;
        catch
            warning('error encountered running config %d',i);
        end
    end

    if saveconfigs
        safeconfigs = configs;
        save('safeconfigs.mat','safeconfigs');
    end
else
    if do_rationalize
        load('safeconfigs.mat');
        configs = safeconfigs;
    else
        load('ratconfigs.mat');
        configs = ratconfigs;
    end
end

%% Rationalize Coefficients
if do_rationalize
    digits(ndigits);
    ratconfigs = configs;
    parfor i=1:length(ratconfigs)
        digits(ndigits);
        config = ratconfigs{i};
        isbtarget = .9999*config.isbn;
        numapprox = config.numapprox;
        results = reducerats(config.cf,config.nfs,config.nds,1e-6,4,numapprox,isbtarget,0);

        results = sortbyfield(results,'ndigits_dep');
        results = sortbyfield(results,'ndigits_free');

        idx = 0;
        isbn = 0;
        while idx < length(results) && isbn < isbtarget
            idx = idx+1;
            rc = results{idx};
            isbn = double(isb_vpa(rc.cf,rc.nfs,rc.nds)/(rc.N+1));
            rc.isbn = isbn;
            ratconfigs{i} = rc;
        end
        if idx > length(results)
            warning('No suitable rational approximation found!');
        end
    end
    if saveconfigs
        save('ratconfigs.mat','ratconfigs');
    end
    configs = ratconfigs;
end

%% Plot
for i=1:length(configs)
    order = configs{i}.p;
    ncores = configs{i}.ncores;
    isbn = double(configs{i}.isbn);
    optpol = double(configs{i}.optpol);
    
    % First subplot
    fig = figure;
    set(fig,'Position',[100 100 1000 600]);
    set(fig,'Color','w');
    ax = subplot(1,2,1);
    
    % Complex xi grid
    xlimits = [-.4 .05];
    ylimits = [-1.05 1.05];
    npoints = 2048;
    [xx,yy] = meshgrid(linspace(xlimits(1),xlimits(2),npoints),linspace(ylimits(1),ylimits(2),npoints));
    xi = complex(xx,yy);

    % Plot
    plot_stability_domain(ax,optpol,xi);
    h = ax; h.FontSize=16;
    h = xlabel('Re{(\xi)}'); h.FontSize=22;
    h = ylabel('Im{((\xi)}'); h.FontSize=22;
    axis(ax,[xlimits ylimits]);

    % Annotation
    pos = get(ax,'Position');
    isbprint = floor(isbn*1e4)/1e4;
    legendtext = sprintf('  %dth Order,   ISB_n = %1.4f',order,isbprint);
    annotation('textbox',[pos(1)+.025*pos(3),pos(2)+.89*pos(4),.068,.065],'String',legendtext,...
                     'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                     'EdgeColor',[0 0 0],'FaceAlpha',1,'FontSize',18);
    
    % Second subplot
    ax = subplot(1,2,2);

    % Complex xi grid
    xlimits = [-.001 .006];
    ylimits = [.15 .85];
    npoints = 2048;
    [xx,yy] = meshgrid(linspace(xlimits(1),xlimits(2),npoints),linspace(ylimits(1),ylimits(2),npoints));
    xi = complex(xx,yy);

    % Plot
    plot_stability_domain(ax,optpol,xi);
    h = ax; h.FontSize=16;
    h = xlabel('Re{(\xi)}'); h.FontSize=18;
    h = ylabel('Im{((\xi)}'); h.FontSize=18;
    axis(ax,[xlimits ylimits]);

    % Export
    paper_export_fig(fig,sprintf(sprintf('stabdomain_p%d_%dcore_rat',order,ncores)));
end

function lam = make_curve(npoints, yoffset, alphaup, alphadn, scale)
    assert(yoffset+alphaup+alphadn <= 1);
    lam = 1i*linspace(0,1,npoints);

    offset = find(imag(lam)>yoffset,1);
    indup = find(imag(lam)>yoffset+alphaup,1);
    inddn = find(imag(lam)>1-alphadn,1);
    nup = indup-offset;
    ndn = length(lam)-inddn;
    nflat = npoints-offset-nup-ndn;
    lam(offset+1:offset+nup) = lam(offset+1:offset+nup)+scale*.5*(1-cos(linspace(0,1,nup)*pi));
    lam(offset+nup+1:end-ndn) = lam(offset+nup+1:end-ndn)+scale*ones(1,nflat);
    lam(end-ndn+1:end) = lam(end-ndn+1:end)+scale*cos(linspace(0,1,ndn)*pi/2);
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
