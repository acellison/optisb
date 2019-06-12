function plot_stability_domain(ax,pol,zgrid,plotaxes)
    sz = size(pol);
    if sz(1) > sz(2)
        pol = pol';
    end
    npoints = 512;
    if nargin >= 3 && ~isempty(zgrid)
        xx = real(zgrid);
        yy = imag(zgrid);
        xi = zgrid;
    else
        [xx,yy] = meshgrid(linspace(-1.05,0.1,npoints),linspace(-1.05,1.05,npoints));
        xi = complex(xx,yy);
    end
    if nargin < 4
        plotaxes = 1;
    end
    fnz = find(pol~=0,1);
    order = length(pol)-fnz;
    pev = polyval(double(pol),order*xi);
    contour(ax,xx,yy,abs(pev),[1 1],'LineColor',[0, 0.4470, 0.7410],'LineWidth',1.8);
    hold on;
    if plotaxes
        l1 = plot(ax,[0 0],[-1.05 1.05],'k--','LineWidth',1.0,'HandleVisibility','off');
        uistack(l1,'bottom');
        l2 = plot(ax,[-1. .1],[0 0],'k--','LineWidth',0.8,'HandleVisibility','off');
        uistack(l2,'bottom');
    end
    grid on; 
    %    plot(ax,real(zers),imag(zers),'.k');
end
