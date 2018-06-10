function h = scatterhist(x,y,varargin)
%SCATTERHIST 2D scatter plot with marginal histograms.
%   SCATTERHIST(X,Y) creates a 2D scatterplot of the data in the vectors X
%   and Y, and puts a univariate histogram on the horizontal and vertical
%   axes of the plot.  X and Y must be the same length.
%
%   SCATTERHIST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values to control how the plot is made. Valid
%   parameters are the following:
%
%       'NBins'    - A scalar or a two-element vector specifying the number
%                    of bins for the X and Y histograms.  The default is to
%                    compute the number of bins using Scott's rule based on
%                    the sample standard deviation.
%
%      'Location'  - A string controlling the location of the marginal
%                    histograms within the figure.  'SouthWest' (the
%                    default) plots the histograms below and to the left of
%                    the scatterplot, 'SouthEast' plots them below and to
%                    the right, 'NorthEast' above and to the right, and
%                    'NorthWest' above and to the left.
%
%      'Direction' - A string controlling the direction of the marginal
%                    histograms in the figure. 'in' (the default) plots the
%                    histograms with bars directed in towards the
%                    scatterplot, 'out' plots the histograms with bars
%                    directed out away from the scatterplot.
%
%      'Group'     - A grouping variable. SCATTERHIST(X,Y,'Group',G)
%                    creates a 2D GSCATTER plot instead of a SCATTER plot,
%                    and the marginal histograms are replaced by grouped
%                    kernel density plots.
%                    
%      'Kernel'    - Indicates if histograms are replaced by grouped kernel
%                    density plots. Choices are:
%        'on'        - This is the default when the 'Group' option is used. 
%        'off'       - Histograms are used. This is the default if no
%                      grouping varialbe is specified, and grouping is
%                      ignored if it is specified.
%        'overlay'   - A kernel density plot of the overall marginal
%                      distribution (grouping is ignored) is overlaid on
%                      the histogram, similar to HISTFIT(x,'kernel').
%
%      'Bandwidth' - The bandwidth of the kernel smoothing window. The
%                    default is empty and SCATTERHIST finds the optimal
%                    bandwidths for estimating normal densities. However,
%                    you may want to choose smaller values to reveal small
%                    features, and it must be a matrix of size 2 by K,
%                    where K is the number of unique groups. The first row
%                    gives the bandwidth of each group in X, and the second
%                    row gives the bandwidth of each group in Y.
%
%      'Legend'    - Indicates if a legend is created when a grouping
%                    variable is given. Choices:
%       'on' or true     - Sets legend visible. This is the default when a
%                          grouping variable is specified.
%       'off' or false   - Sets legend invisible. This is the default when
%                          no grouping variable is specified
%
%   The following parameters have no effect on the histogram when 'Kernel'
%   is set to either 'off' or 'Overlay'.
%      'LineStyle' - A cell vector of valid line style strings that specify
%                    the 'LineStyle' property of each kernel density line.
%                    For example, {'-',':','-.'}. See documentation of PLOT
%                    for a list of valid line style strings.
%                   
%      'LineWidth' - A scalar or a vector of numbers that specify the
%                    'LineWidth' property of each kernel density line. For
%                    example, [1 2 3].
%
%      'Color'     - A string or a matrix of RGB values that specifies the
%                    marker color of each group in the gscatter plot. For
%                    example, 'r' or {'red'} if not grouped, or 'rbk' if
%                    grouped. See the ColorSpec documentation for a table
%                    that lists the predefined colors and their RGB
%                    equivalents.
% 
%      'Marker'    - A string that specifies the markers used for each
%                    group in the gscatter plot. See the documentation of
%                    PLOT for a list of valid markers. The default is 'o',
%                    a circle.
%
%      'MarkerSize'- A numerical vector that specifies the size of the each
%                    marker. The default value is 6.
%
%   NOTE:SCATTERHIST cycles through the available values of the parameters:
%   'LineWidth', 'LineStyle', 'Color', 'Marker' and 'MarkerSize' when the
%   total number of groups is greater than the number of specified values.
%
%
%   H = SCATTERHIST(...) returns a vector of three axes handles for the
%   scatterplot, the plot along the horizontal axis, and the plot along the
%   vertical axis, respectively.
%
%   SCATTERHIST(X,Y,NBINS) is supported for backwards compatibility.
%
%   Example:
%      Independent normal and lognormal random samples
%         x = randn(1000,1);
%         y = exp(.5*randn(1000,1));
%         scatterhist(x,y)
%      Marginal uniform samples that are not independent
%         u = copularnd('Gaussian',.8,1000);
%         scatterhist(u(:,1),u(:,2))
%      Mixed discrete and continuous data
%         load('carsmall');
%         scatterhist(Weight,Cylinders,'NBins',[10 3])
%
%   See also gscatter, plot, hist, histfit, grpstats, grp2idx, ksdensity

%   Copyright 2006-2012 The MathWorks, Inc.


narginchk(2, Inf);

if ~isvector(x) || ~isnumeric(x) || ~isvector(y) || ~isnumeric(y)
    error(message('stats:scatterhist:BadXY'));
end
if numel(x)~=numel(y)
    error(message('stats:scatterhist:BadXY'));
end
x = x(:);
y = y(:);

if nargin == 3
    h = scatterhist(x,y,'nbins',varargin{1});
    return;
else
    % Parse name/value pairs
    pnames = {'nbins','location','direction','group','kernel',...
        'color', 'marker', 'markersize','linestyle','linewidth',...
        'legend', 'bandwidth', 'Xfunc', 'Yfunc', 'Xlim', 'Ylim'};
    
    dflts =  {  [],  'SouthWest',  'in',   [],   'off' ...
        [],   'o'         6        {'-',':','-.','--'}, 1, ...
        'on',        [],           []        [], [], []};
    
    [nbins,location,direction,group,kernel,clr,marker, ...
        siz, lStyle, lw,doLegend, ww, Xfunc, Yfunc,Xlim,Ylim,setFlag]...
        = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    % Validate values
    if setFlag.nbins &&...
            (numel(nbins)>2 ||~internal.stats.isIntegerVals(nbins))
        error(message('stats:scatterhist:BadBins'));
    end   
                                          
    kernelKWDs = {'on','off','overlay'};
    kernel = internal.stats.getParamVal(kernel,kernelKWDs,'''Kernel''');
    
    if ~(isnumeric(siz)&& all(siz>0))
        error(message('stats:scatterhist:BadNumericLineProperty',...
            '''MarkerSize'''));
    end
    
    if ~iscellstr(lStyle)
        error(message('stats:scatterhist:BadLineStyle'));
    end
    if ~(isnumeric(lw)&& all(lw>0))
        error(message('stats:scatterhist:BadNumericLineProperty',...
            '''LineWidth'''));
    end
    
    doLegend = internal.stats.parseOnOff(doLegend,'''Legend''');
end

if ~isempty(group)    
    [gInt,gn,~] = grp2idx(group);
    if numel(gInt) ~= numel(x)
        error(message('stats:scatterhist:BadGroup'));
    end
    k = numel(gn); % number of unique groups
    
    if ~setFlag.kernel
        kernel = 'on'; % default 'on' if there is grouping variable
    end
    
    if ~setFlag.legend
        doLegend = true;
    end
    
    [marker,siz,lStyle,lw] = ...
        internal.stats.cycleLineProperties(k,marker,siz,lStyle,lw);
else
    k = 1;
end

if isempty(clr)
    clr = lines(k);
else
    if isnumeric(clr) && isrow(clr)
        if length(clr) == 3
        % if it is a row vector, replicate it to a 2-row matrix and pass to
        % internal.stats.cycleLineProperties
            clr = [clr;clr];
        else
            error(message('stats:internal:colorStringToRGB:ValueMustBe3ElementVector'));
        end
    end
    clr = internal.stats.cycleLineProperties(k,clr);
    clr = internal.stats.colorStringToRGB(clr);
end

if isempty(nbins)
        % By default use the bins given by Scott's rule
        [xctrs,yctrs] = defaultBins(x(~isnan(x)),y(~isnan(y)));
    elseif isscalar(nbins)
        xctrs = nbins;  % specified number of bins, same for x and y
        yctrs = nbins;
    else
        xctrs = nbins(1); % specified number of bins, x and y different
        yctrs = nbins(2);
end

% Set the bandwidth of kernel smoothing window
if isempty(ww) 
    ww1 = [];
    ww2 = [];
elseif isnumeric(ww) && isequal(size(ww),[2,k]) && (all(~isnan(ww(:))))
    ww1 = ww(1,:);
    ww2 = ww(2,:);
else
    error(message('stats:scatterhist:BadBandwidth',k));
end

% Setup marginal distribution plotting function
if isempty(Xfunc) || isempty(Yfunc)
    switch kernel
        case 'on'
            if isempty(group)
                group = ones(1,numel(x));
                doLegend = false;
            end
            Xfunc = @()internal.stats.plotGroupedKSDensity(x,group,...
                'Color',clr,'LineWidth',lw,'LineStyle',lStyle,'Width',ww1);
            Yfunc = @()internal.stats.plotGroupedKSDensity(y,group,...
                'Color',clr,'LineWidth',lw,'LineStyle',lStyle,'Width',ww2);
            
        case 'overlay'
            if ~setFlag.nbins
                xctrs = numel(xctrs);
                yctrs = numel(yctrs);
            end
            Xfunc = @()histfit(x,xctrs,'kernel');
            Yfunc = @()histfit(y,yctrs,'kernel');
          
        case 'off'
            Xfunc = @()hist(x,xctrs);
            Yfunc = @()hist(y,yctrs);
            
    end
end

% Set up positions for the plots
switch lower(direction)
    case 'in'
        inoutSign = 1;
    case 'out'
        inoutSign = -1;
    otherwise
        error(message('stats:scatterhist:BadDirection'));
end

switch lower(location)
    case {'ne' 'northeast'}
        scatterLoc = 3;
        scatterPosn = [.1 .1 .55 .55];
        scatterXAxisLoc = 'top'; scatterYAxisLoc = 'right';
        histXLoc = 1; histYLoc = 4;
        histXSign = -inoutSign;
        histYSign = -inoutSign;
    case {'se' 'southeast'}
        scatterLoc = 1;
        scatterPosn = [.1 .35 .55 .55];
        scatterXAxisLoc = 'bottom'; scatterYAxisLoc = 'right';
        histXLoc = 3; histYLoc = 2;
        histXSign = inoutSign;
        histYSign = -inoutSign;
    case {'sw' 'southwest'}
        scatterLoc = 2;
        scatterPosn = [.35 .35 .55 .55];
        scatterXAxisLoc = 'bottom'; scatterYAxisLoc = 'left';
        histXLoc = 4; histYLoc = 1;
        histXSign = inoutSign;
        histYSign = inoutSign;
    case {'nw' 'northwest'}
        scatterLoc = 4;
        scatterPosn = [.35 .1 .55 .55];
        scatterXAxisLoc = 'top'; scatterYAxisLoc = 'left';
        histXLoc = 2; histYLoc = 3;
        histXSign = -inoutSign;
        histYSign = inoutSign;
    otherwise
        error(message('stats:scatterhist:BadLocation'));
end

% Create empty subplots
if ~matlab.graphics.internal.isGraphicsVersion1() 
    clf;
end
hScatter=subplot(2,2,scatterLoc);
hYHist=subplot(2,2,histYLoc);
hXHist=subplot(2,2,histXLoc);

% Use the function handles Xfunc and Yfunc to plot the marginal
% distribution
%colormap([0.8 0.8 1]); % set the default face color of the bars

% 'xfunc' and 'yfunc' are NOT external parameters. The function handles can
% reference functions of two kinds of signatures. The second signature may
% change in the future.
set(gcf,'CurrentAxes',hXHist);
if nargin(Xfunc) ==0
    Xfunc();
else 
    Xfunc(x,group,'Color',clr,'LineStyle',lStyle);
end

set(gcf,'CurrentAxes',hYHist);
if nargin(Yfunc)==0
    Yfunc();
else
    Yfunc(y,group,'Color',clr,'LineStyle',lStyle);
end

view(hYHist,270,90); % Rotate the yhist plot
if histXSign<0
    set(hXHist,'YDir','reverse');    
end
if histYSign>0
    set(hYHist,'YDir','reverse');
end

% Scatter plot by gscatter 
set(gcf,'CurrentAxes',hScatter);
if isempty(inputname(1))
    xname = 'x';
else
    xname = inputname(1);
end
if isempty(inputname(2))
    yname = 'y';
else
    yname = inputname(2);
end
%gscatter(x,y,group,clr,marker,siz,'on',xname,yname);
nbins=100;
allD=[x, y];
if length(Xlim)>0
    x_bin_size=(Xlim(2)-Xlim(1))/nbins;
    x_ctrs=linspace(Xlim(1)+x_bin_size/2,Xlim(2)-x_bin_size/2,nbins);
else
    x_bin_size=(max(x)-min(x))/nbins;
    x_ctrs=linspace(min(x)+x_bin_size/2,max(x)-x_bin_size/2,nbins);
end
if length(Ylim)>0
    y_bin_size=(Ylim(2)-Ylim(1))/nbins;
    y_ctrs=linspace(Ylim(1)+y_bin_size/2,Ylim(2)-y_bin_size/2,nbins);
else
    y_bin_size=(max(y)-min(y))/nbins;
    y_ctrs=linspace(min(y)+y_bin_size/2,max(y)-y_bin_size/2,nbins);
end
allN=[];
for s_idx=1:length(unique(group))
    idx=find(group==s_idx);
    [N,c]=hist3(allD(idx,:),'Ctrs',{x_ctrs y_ctrs});
    x_diff=diff(x_ctrs);
    y_diff=diff(y_ctrs);
    x = repmat([x_diff, x_diff(end)], length(x_ctrs),1)';
    y = repmat([y_diff, y_diff(end)], length(y_ctrs),1);
    % volume of the histogram
    V_tot  = sum(sum(x.*y.*N));
    allN(s_idx,:,:)=N./V_tot.*100.0;
end
imagesc(x_ctrs,y_ctrs,squeeze(mean(allN))','Parent',hScatter);
set(hScatter,'ydir','normal');
colorbar('WestOutside');

% Create Legend and set visibility accordingly
hLegend = findobj(gcf,'Tag','legend');
hLegendToggle = findall(gcf,'Tag','Annotation.InsertLegend');
if ~isempty(hLegend)&& ~doLegend
    set(hLegend,'Visible','off');
    set(hLegendToggle,'State','off');
end
localSetLegend();

% Change legend toggle button call back functions
set(hLegendToggle,'ClickedCallback',@legendtogglecallback);

% adjust axes settings
set(hScatter,'Box','on');
axis(hScatter,'tight');
axis(hXHist,'off','tight');
axis(hYHist,'off','tight');

% Make sure that the axes do not overlay or overlap on the plots so that
% the edges of the xhist and yhist are not being cropped by the HG renderer
bestLimFunc = @(x,y)[min(x(1),y(1)), max(x(2),y(2))];
%XhistXlim = bestLimFunc(get(hXHist,'Xlim'),get(hScatter,'Xlim'));
%YhistXlim = bestLimFunc(get(hYHist,'Xlim'),get(hScatter,'Ylim'));
%set(hScatter,'Xlim', XhistXlim + [-0.01, 0.01]*range(XhistXlim));
%set(hScatter,'Ylim', YhistXlim + [-0.01, 0.01]*range(YhistXlim));
set(hXHist,'XLim',get(hScatter,'XLim')); 
set(hYHist,'XLim',get(hScatter,'YLim')); 

XhistYlim = get(hXHist,'YLim');
YhistYlim = get(hYHist,'YLim');
set(hXHist,'YLim',XhistYlim + [-0.01, 0.01]*range(XhistYlim)); 
set(hYHist,'YLim',YhistYlim + [-0.01, 0.01]*range(YhistYlim));

% Invisible texts to help calculate positions
txt1 = text(0,0,'','Visible','off','HandleVisibility','off');
txt2 = text(1,1,'','Visible','off','HandleVisibility','off');

set(hScatter,'Position',scatterPosn, 'XAxisLocation',scatterXAxisLoc, ...
             'YAxisLocation',scatterYAxisLoc, 'Tag','scatter');
set(hXHist,'Tag','xhist');
set(hYHist,'Tag','yhist');

% Sync axis
myscatterhistXYLimCallback();

% Add listeners to resize or relimit histograms when the scatterplot changes
flagUpdate = true;
%addmyscatterhistlisteners(hScatter, @scatterhistPositionCallback,...
%       @scatterhistXYLimCallback, @markedCleanCallback);

% Disable unrelated toolbar ui buttons
internal.stats.removeFigToolbarButton(gcf,{'Exploration.Brushing',...
                                           'Annotation.InsertColorbar',...
                                           'Exploration.Rotate',...
                                           'DataManager.Linking'});
% Disable pan on xhist and yhist plots
hpan = pan;
hpan.setAllowAxesPan([hXHist,hYHist],false);
% Disable zoom on xhist and yhist plots
hzoom = zoom;
hzoom.setAllowAxesZoom([hXHist,hYHist],false);

% Set current axis to the scatter axis
scatterhistPositionCallback();
set(gcf,'CurrentAxes',hScatter);
set(gcf,'NextPlot','replace');

if nargout>0
    h  = [hScatter hXHist hYHist];
end

%--------------------------------------------------------------------------
% Nested subfunctions
    function markedCleanCallback(~,~)
        flagUpdate = true;
    end

    function scatterhistPositionCallback(~,~)
    % Position Callback Function
    
    if ~flagUpdate
        return;
    end
    
    posn = getrealposition(hScatter,txt1,txt2);
    if ~matlab.graphics.internal.isGraphicsVersion1()
        % drawnow ensures that OuterPosition is updated
        drawnow('expose');
    end
    oposn = get(hScatter,'OuterPosition');
    posnLegend = get(hLegend,'Position');

    switch lower(location)
        case {'sw' 'southwest'}
            % vertically: margin, histogram, margin/4, scatterplot, margin
            vmargin = min(max(1 - oposn(2) - oposn(4), 0), oposn(2));
            posnHistX = [posn(1) vmargin posn(3) oposn(2)-1.25*vmargin];
            % horizontally: margin, histogram, margin/4, scatterplot, margin
            hmargin = min(max(1 - oposn(1) - oposn(3), 0), oposn(1));
            posnHistY = [hmargin posn(2) oposn(1)-1.25*hmargin posn(4)];
            % legend position
            if ~isempty(posnLegend)
                posnLegend = [max(oposn(1)-posnLegend(3),0.1*hmargin),...
                    max(oposn(2)-posnLegend(4),0.1*vmargin),...
                    posnLegend(3),posnLegend(4)];
            end
        case {'ne' 'northeast'}
            % vertically: margin, scatterplot, margin/4, histogram, margin
            vmargin = max(oposn(2), 0);
            posnHistX = [posn(1) oposn(2)+oposn(4)+.25*vmargin ...
                posn(3) 1-oposn(2)-oposn(4)-1.25*vmargin];
            % horizontally: margin, scatterplot, margin/4, histogram, margin
            hmargin = max(oposn(1), 0);
            posnHistY = [oposn(1)+oposn(3)+.25*hmargin    posn(2)...
                1-oposn(1)-oposn(3)-1.25*hmargin  posn(4)];
            % legend position
            if ~isempty(posnLegend)
                posnLegend = [...
                    min(oposn(1)+oposn(3),1-posnLegend(3)-0.1*hmargin),...
                    min(oposn(2)+oposn(4),1-posnLegend(4)-0.1*hmargin),...
                    posnLegend(3),...
                    posnLegend(4)];
            end
        case {'se' 'southeast'}
            % vertically: margin, histogram, margin/4, scatterplot, margin
            vmargin = max(1 - oposn(2) - oposn(4), 0);
            posnHistX = [posn(1) vmargin posn(3) oposn(2)-1.25*vmargin];
            % horizontally: margin, scatterplot, margin/4, histogram, margin
            hmargin = max(oposn(1), 0);
            posnHistY = [oposn(1)+oposn(3)+.25*hmargin     posn(2)...
                1-oposn(1)-oposn(3)-1.25*hmargin  posn(4)];
            % legend position
            if ~isempty(posnLegend)
                posnLegend = [...
                    min(oposn(1)+oposn(3), 1-posnLegend(3)-0.1*hmargin),...
                    max(oposn(2)-posnLegend(4),0.1*vmargin),...
                    posnLegend(3),...
                    posnLegend(4)];
            end
        case {'nw' 'northwest'}
            % vertically: margin, scatterplot, margin/4, histogram, margin
            vmargin = max(oposn(2), 0);
            posnHistX = [posn(1) oposn(2)+oposn(4)+.25*vmargin...
                posn(3) 1-oposn(2)-oposn(4)-1.25*vmargin];
            % horizontally: margin, histogram, margin/4, scatterplot, margin
            hmargin = max(1 - oposn(1) - oposn(3), 0);
            posnHistY = [hmargin posn(2) oposn(1)-1.25*hmargin posn(4)];
            % legend position
            if ~isempty(posnLegend)
                posnLegend = [max(oposn(1)-posnLegend(3),0.1*hmargin),...
                    min(oposn(2)+oposn(4),1-posnLegend(4)-0.1*hmargin),...
                    posnLegend(3),...
                    posnLegend(4)];
            end
    end
    posnHistX = max(posnHistX,[0 0 .05 .05]);
    posnHistY = max(posnHistY,[0 0 .05 .05]);

    set(hXHist,'Position',posnHistX);
    set(hYHist,'Position',posnHistY);
    if ~isempty(posnLegend)
        set(hLegend, 'Position',posnLegend);
    end

    myscatterhistXYLimCallback();
    
    if ~matlab.graphics.internal.isGraphicsVersion1()
        flagUpdate = false;
    end
    end % scatterhistPositionCallback function

% -----------------------------
    function myscatterhistXYLimCallback(~,~)
    %callback function to sync the x axes of xhist and yhist with the axes
    %of scatter plot.

    set(hXHist,'XLim',get(hScatter,'XLim'));
    set(hYHist,'XLim',get(hScatter,'YLim'));
    end


% ------------------------------
    function localSetLegend(~,~)
    if isempty(lStyle) % If no line object, no line added to legends
        return;
    elseif numel(findobj(hXHist,'Type','line','-depth',1)) <3
        % there should be more than 2 lines
        return;
    end

    hLegLines = findobj(hLegend,'Type','line');
    for i = 0:2:(numel(hLegLines)-2)
        set(hLegLines(end-i),'LineStyle',lStyle{i/2+1},...
            'LineWidth',lw(i/2+1));
    end
    end

% ------------------------------
    function legendtogglecallback(~,~)
    if strcmp(get(hLegend,'Visible'),'off')
        set(hLegend,'Visible','on')
    else
        set(hLegend,'Visible','off')
    end
    end

end % SCATTERHIST main function


% SUBFUNCTIONS:
% -----------------------------
function p = getrealposition(a,txt1,txt2)
p = get(a,'Position');

% For non-warped axes (as in "axis square"), recalculate another way
if isequal(get(a,'WarpToFill'),'off')
    pctr = p([1 2]) + 0.5 * p([3 4]);
    xl = get(a,'XLim');
    yl = get(a,'YLim');
    
    % Use text to get coordinate (in points) of southwest corner
    set(txt1,'Units','data');
    set(txt1,'Position',[xl(1) yl(1)]);
    set(txt1,'Units','pixels');
    pSW = get(txt1,'Position');
    
    % Same for northeast corner
    set(txt2,'Units','data');
    set(txt2,'Position',[xl(2) yl(2)]);
    set(txt2,'Units','pixels');
    pNE = get(txt2,'Position');
    
    % Re-create position
    % Use min/max/abs in case one or more directions are reversed
    p = [min(pSW(1),pNE(1)), ...
         max(pSW(2),pNE(2)), ...
         abs(pNE(1)-pSW(1)), ...
         abs(pNE(2)-pSW(2))];
    p = hgconvertunits(ancestor(a,'figure'),p, ...
             'pixels','normalized',ancestor(a,'figure'));
    
    % Position to center
    p = [pctr(1)-p(3)/2, pctr(2)-p(4)/2, p(3), p(4)];
end
end % getrealposition function


function [xctrs,yctrs] = defaultBins(x,y)
% By default use the bins given by Scott's rule
xctrs = internal.stats.histbins(x); % returns bin ctrs
yctrs = internal.stats.histbins(y); % returns bin ctrs
if length(xctrs)<2
    xctrs = 1; % for constant data, use one bin
end
if length(yctrs)<2
    yctrs = 1; % for constant data, use one bin
end
end % defaultBins function
