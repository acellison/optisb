function paper_export_fig(ax,filename,linewidth)
    if nargin<3
        linewidth = 1.4;
    end
    
    basepath = 'paper_figures/';
    epspath = [basepath 'eps/'];
    graypath = [epspath 'gray/'];
    if ~exist(basepath,'dir')
        mkdir(basepath)
    end    
    if ~exist(epspath,'dir')
        mkdir(epspath)
    end    
    if ~exist(graypath,'dir')
        mkdir(graypath)
    end    

    set(gcf,'Renderer','painters');
    saveas(gcf,[epspath filename],'epsc');

    if strcmp(ax.Type,'figure')
        ch = allchild(ax);
        for i=1:length(ch)
            if strcmp(ch(i).Type,'axes')
                set(ch(i),'color',[1,1,1]);
                tobw(ch(i),linewidth);
            end
        end
    else
        tobw(ax,linewidth);
    end
    
    set(gcf,'Renderer','painters');
    saveas(gcf,[graypath filename],'epsc');
end

function tobw(ax,linewidth)
    ch = allchild(ax);
    for i=1:length(ch)
        graphic = ch(i);

        check_and_set(graphic,'LineColor','k');
        check_and_set(graphic,'Color','k');

        % Axes lines should be less thick
        lw = linewidth;
        if strcmp(graphic.Type,'line')
            if all(graphic.XData==0)
                lw = 1.0;
            elseif all(graphic.YData==0)
                lw = 0.8;
            end
            if length(graphic.XData)==1
                lw = 0.7;
            end
        end
        check_and_set(graphic,'LineWidth',lw);
    end
end

function check_and_set(h,field,value)
    if isprop(h,field) || isfield(h,field)
        set(h,field,value);
    end
end