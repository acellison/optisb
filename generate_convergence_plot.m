clear; close all;

%% User Settings
do_generate = 1;

% Can be one of 'double' or 'mp'.  If set to 'mp', will use the
% Advanpix Multiprecision Computing toolbox for all PDE data.  This yields
% the nicest convergence plots since we can get way down below 10^{-16}.
% Since Matlab's fft routine doesn't support 'vpa' this data type is 
% unsupported for now.
global hp_class;
hp_class = 'mp';

%% Setting-dependent variables
use_mp = strcmp(hp_class,'mp');
dataname = sprintf('convect_configs_%s.mat',hp_class);

%% Generate
if do_generate
    ndigits = 100;
    setdigits(ndigits)

    %% Set up the solver
    npoints = 2;                                 % Points on the grid
    kmax = 1;                                     % Max Grid size
    c = 1;                                        % Wave speed
    ugen = @(xi)0.5*(1-cos(2*high_precision('pi')*xi));
    spectral = 1;                                 % Use spectral derivatives
    hscales = high_precision(2.^(-4:-1:-8));
    gen = @(t,npoints)generator(ugen,t,c,npoints);

    %% Set up the schemes
    configs = cell(1,1); i=0;

    i=i+1; configs{i} = struct('p',8,'isbn',.7675, ...
        'nds',[2 4 6 10],'nfs',[8 12 14 16 18 20 22], ...
        'cfs', [sym('269/95360') sym('13805/611712') sym('83/1314') sym('3607/16544') ...
                sym('2196/619') sym('-35962/2395') sym('7352/605')]);
    i=i+1; configs{i} = struct('p',12,'isbn',.7116, ...
        'nds',[2 8 10 16 24 26],'nfs',[4 6 12 14 18 20 22 28 30], ...
        'cfs', [sym('12985/994150711296') sym('3295/1296039936') sym('2521/8515584')  sym('1349/1959936') ...
                sym('11223/1226368') sym('5711/69600') sym('6007/2478') sym('-1338112/5553') sym('50764/475')]);

    for i=1:length(configs)
        nfs = configs{i}.nfs(:);
        nds = configs{i}.nds(:);
        cfs = configs{i}.cfs(:);
        Vdep = rot90(vander(1./sym(nds).^2));
        Vfree = rot90(vandern(1./sym(nfs).^2,length(nds)));
        cds = Vdep\([1;zeros(length(nds)-1,1)]-Vfree*cfs);
        configs{i}.nfs = nfs;
        configs{i}.nds = nds;

        if use_mp
            configs{i}.cfs = sym2mp(cfs);
            configs{i}.cds = sym2mp(cds);
        else
            configs{i}.cfs = cfs;
            configs{i}.cds = cds;
        end
        configs{i}.N = max(max(nfs),max(nds));
    end

    %% Set the time step size
    normfn = @(x,y)max(abs(x-y));
    tn = high_precision(npoints);  % one revolution
    t0 = high_precision(0);

    %% Run the solvers
    for i=1:length(configs)
        setdigits(ndigits);

        order = configs{i}.p;
        maxthread = configs{i}.N;
        isbn = configs{i}.isbn;

        ns = [configs{i}.nds;configs{i}.nfs];
        cs = [configs{i}.cds;configs{i}.cfs];
        cs = high_precision(cs);

        hmax = .99*kmax*(maxthread+1)*isbn/c;
        if spectral, hmax = hmax/pi; end

        %% Revolve and compute error on the final revolution
        num_measurements = length(hscales);
        errors = zeros(num_measurements,1,class(high_precision(0)));
        for j=1:num_measurements
            hcur = hscales(j)*hmax;
            kcur = hscales(j)*kmax;
            ncur = ceil((tn-t0)/hcur);
            npointscur = round(npoints/kcur);

            u0 = gen(t0,npointscur);
            wave = @(t,u)onewaywave(t,u,kcur,c,spectral);
            un = rexstep(@gbs,wave,hcur,t0,u0,ns,cs,ncur);

            utrue = gen(t0+ncur*hcur/kcur,npointscur);
            errors(j) = normfn(un,utrue);
        end

        configs{i}.errors = errors;
        configs{i}.hs = hmax*hscales;
        configs{i}.hscaled = configs{i}.hs/(maxthread+1);
    end

    %% RK4
    rkhmax = .99*kmax*sqrt(8)/c;
    if spectral, rkhmax = rkhmax/pi; end

    % Revolve and compute error on the final revolution
    num_measurements = length(hscales);
    rkerrors = zeros(num_measurements,1,class(high_precision(0)));
    for j=1:num_measurements
        hcur = hscales(j)*rkhmax;
        kcur = hscales(j)*kmax;
        ncur = ceil((tn-t0)/hcur);
        npointscur = round(npoints/kcur);

        u0 = gen(t0,npointscur);
        wave = @(t,u)onewaywave(t,u,kcur,c,spectral);
        un = rk4(wave,hcur,t0,u0,ncur);

        utrue = gen(t0+ncur*hcur/kcur,npointscur);
        rkerrors(j) = normfn(un,utrue);
    end
    configs{3} = struct('kind','rk','p',4,'errors',rkerrors,'hs',rkhmax*hscales, ...
                        'hscaled',rkhmax*hscales/4);

	if use_mp
        for i=1:2
            configs{i}.cfs = mp2sym(configs{i}.cfs);
            configs{i}.cds = mp2sym(configs{i}.cds);
        end
        for i=1:length(configs)
            configs{i}.errors = mp2sym(configs{i}.errors);
            configs{i}.hs = mp2sym(configs{i}.hs);
            configs{i}.hscaled = mp2sym(configs{i}.hscaled);
        end
    end
    save(dataname,'configs');
else
    load(dataname);
end


%% Plot error
figure;
set(gcf,'Color','w');
set(gcf,'renderer','Painters');
colors = linecolors();
linewidth = 1.4;
loglog(configs{end}.hscaled, configs{end}.errors,'-*','LineWidth',linewidth,'Color',colors(1,:));
hold on;
markers = {'-v','-s','-^','-o'};
for i=1:length(configs)-1
    loglog(configs{i}.hscaled, configs{i}.errors,markers{i},'LineWidth',linewidth,'Color',colors(i+1,:));
end
linewidth = 1.2;
plot([2e-3 6e-3],.7e-9*[1 (3)^4],'--k','LineWidth',linewidth);
plot([2e-3 6e-3],.7e-16*[1 (3)^8],'--k','LineWidth',linewidth);
plot([2e-3 6e-3],.7e-24*[1 (3)^12],'--k','LineWidth',linewidth);
grid on;
axis([6e-4 2e-2 1e-29 1e-3]);
h = get(gca,'xaxis'); h.FontSize = 10;
h = get(gca,'yaxis'); h.FontSize = 10;
h = xlabel('Time Step Normalized By Number of Function Evaluations'); h.FontSize=12;
h = ylabel('L_{\infty} Error'); h.FontSize=12;

annotx = .641;
yoff = -.006;
annotation('textbox',[annotx .772+yoff .05 .05],'interpreter','latex', ...
                 'String','$$\mathcal{O}\big({h}^{4}\big)$$', ...
                 'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                 'EdgeColor','none','FaceAlpha',0,'FontSize',12)
annotation('textbox',[annotx .613+yoff .05 .05],'interpreter','latex', ...
                 'String','$$\mathcal{O}\big({h}^{8}\big)$$', ...
                 'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                 'EdgeColor','none','FaceAlpha',0,'FontSize',12)
annotation('textbox',[annotx .423+yoff .05 .05],'interpreter','latex', ...
                 'String','$$\mathcal{O}\big({h}^{12}\big)$$', ...
                 'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                 'EdgeColor','none','FaceAlpha',0,'FontSize',12)

annotx = .24;
annotation('textbox',[annotx .755 .05 .05], ...
                 'String','RK_4', ...
                 'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                 'EdgeColor','none','FaceAlpha',0,'FontSize',12)
annotation('textbox',[annotx .529 .05 .05], ...
                 'String','GBS_{8, 6}', ...
                 'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                 'EdgeColor','none','FaceAlpha',0,'FontSize',12)
annotation('textbox',[annotx .274 .05 .05], ...
                 'String','GBS_{12, 8}', ...
                 'FitBoxToText','on','BackgroundColor',[1 1 1], ...
                 'EdgeColor','none','FaceAlpha',0,'FontSize',12)

paper_export_fig(gca,'convect_error',1);


%% du/dt = -c*du/dx
function ut = onewaywave(~,u,k,c,spectral)
    ux = derivative(u,k,spectral);
    ut = -c*ux;
end

%% Initial data propagator along characteristics x-ct
function u0 = generator(fn,ts,c,npoints)
    u0 = zeros(length(ts),npoints,class(ts));
    for i=1:length(ts)
        u0(i,:) = fn(((0:npoints-1)-c*ts(i))/npoints);
    end
end

%% Derivative routines
function du = derivative(u,k,spectral)
    if spectral
        du = spectral_derivative(u,k);
    else
        du = central_difference(u,k);
    end
end

function du = central_difference(u,k)
    ucirc = [u(end) u u(1)];
    du = (ucirc(3:end)-ucirc(1:end-2))/(2*k);
end

function du = spectral_derivative(u,k)
% function du = spectral_derivative(u,k)
% Compute the p'th derivative of u with grid step size k
    n  = length(u);
    kn = 2.*high_precision('pi')/n/k;
    ik = 1i*kn*[0:n/2 -(n/2)+1:-1];
    du = real(ifft(ik.*fft(u)));
end

%% High precision conversions
function hp = high_precision(x)
    global hp_class;
    if strcmp(hp_class,'double')
        if strcmp(x,'pi')
            hp = pi;
        else
            hp = double(x);
        end
    elseif strcmp(hp_class,'sym')
        hp = vpa(x);
    elseif strcmp(hp_class,'mp')
        hp = mp(x);
    end
end

function setdigits(ndigits)
    global hp_class;
    if strcmp(hp_class,'sym')
        digits(ndigits);
    elseif strcmp(hp_class,'mp')
        mp.Digits(ndigits);
    end
end
