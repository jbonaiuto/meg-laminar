function metric_data=plot_surface_metric(subj_info, surface_file, metric_file, view, varargin)

% Parse inputs
defaults = struct('output_file','','output_format','png','ax',0,'limits',[],'threshold',0,'mask',[],'plot',true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Create figure and axis if not specified
if params.plot || length(params.output_file)>0
    if params.ax==0
        fig=figure('Renderer','OpenGL', 'Color',[1 1 1]);
        params.ax = axes('Visible', 'off', 'Parent', fig, 'Projection', 'perspective', 'PlotBoxAspectRatio', [1.268 1 1.129], 'DataAspectRatio', [1 1 1], 'CameraViewAngle', 6.028, 'CameraUpVector', subj_info.camera_up_vector, 'CameraPosition', subj_info.camera_position);
    % Otherwise set axis properties
    else
        set(params.ax,'Visible','off');
        set(params.ax,'Projection','perspective');
        set(params.ax,'PlotBoxAspectRatio',[1.268 1 1.129]);
        set(params.ax,'DataAspectRatio',[1 1 1]);
        set(params.ax,'CameraViewAngle',6.028);
        set(params.ax,'CameraUpVector',subj_info.camera_up_vector(view));
        set(params.ax,'CameraPosition',subj_info.camera_position(view));
    end

    % Create light
    light('Parent',params.ax,'Style','local','Position',[749 868.1 1263]);

    % Create light
    light('Parent',params.ax,'Style','local','Position',[-1611 -265.3 288.5]);

    % Create light
    light('Visible','off','Parent',params.ax,'Style','local','Position',[776.2 899.5 1309]);
end

% Load surface and metric data
surface=gifti(surface_file);
metric=gifti(metric_file);
metric_data=metric.cdata(:);

if params.plot || length(params.output_file)>0
    % Plot surface
    hp = patch('vertices', surface.vertices, 'faces', surface.faces, 'EdgeColor', 'none', 'Parent', params.ax, 'FaceColor','interp', 'linestyle','none');
    set(hp,'SpecularStrength',0.0);
    set(hp,'AmbientStrength',0.8);

    % Create color maps
    % Control points for positive color map
    pos_control_pts=zeros(3,5);
    pos_control_pts(1,:)=[340 225 90 10 10];
    pos_control_pts(2,:)=[340 340 340 270 10];
    pos_control_pts(3,:)=[340 340 340 340 340];

    % Control points for negative color map
    neg_control_pts=zeros(3,5);
    neg_control_pts(1,:)=[340 264 264 180 340];
    neg_control_pts(2,:)=[10 92 92 264 340];
    neg_control_pts(3,:)=[10 264 264 92 264];

    % If metric value limits are not specified
    if length(params.limits)==0

        % Split metric data into positive and negative values
        pos_vals=find(metric_data>0);
        neg_vals=find(metric_data<0);

        % Compute percentiles
        pos_percentiles=tiedrank(metric_data(pos_vals))/length(pos_vals);
        neg_percentiles=tiedrank(metric_data(neg_vals))/length(neg_vals);

        % Min and max are 2% and 98% of neg and pos values
        min_clipped_val=min(metric_data(neg_vals(find(neg_percentiles>=.02))));
        max_clipped_val=max(metric_data(pos_vals(find(pos_percentiles<=.98))));
        params.limits=[min_clipped_val max_clipped_val];
    end

    % Number of colors overall
    num_colors=255;
    % Compute ratio of positive to negative values
    ratio=abs(params.limits(1))/(abs(params.limits(1))+params.limits(2));
    % Compute number of postive and negative colots
    neglen=round(num_colors*ratio);
    poslen=num_colors-neglen;

    % Generate positive and negative color maps
    C_pos=generate_color_map(poslen, pos_control_pts);
    C_neg=generate_color_map(neglen, neg_control_pts);

    % Set color map
    colormap([C_neg;C_pos]/255);
    % Show colorbar and freeze
    cb=colorbar();
    set(params.ax, 'Clim', params.limits);
    cbfreeze;

    % Generate display version of metric data
    display_metric_data=metric_data;
    % Clip range at the limits
    display_metric_data(find(display_metric_data<params.limits(1)))=params.limits(1);
    display_metric_data(find(display_metric_data>params.limits(2)))=params.limits(2);

    % Calculate color number that each intensity maps to
    numcol = size(colormap,1);
    min_inten = min(display_metric_data);
    cfrac = (display_metric_data - min_inten) ./ (max(display_metric_data) - min_inten);
    display_metric_data = 1 + floor(numcol * (cfrac * (1-2*eps)));
    % Set mapping to direct - intensity values are index of colors in map
    set(hp,'CDataMapping', 'direct');

    % Add extra color (gray to the end of the map)
    colormap([colormap; .5 .5 .5]);
    % Get new size of color map
    numcol = size(colormap,1);
end

% Apply mask and/or threshold
if length(params.mask)>0 || abs(params.threshold)>0

    % Apply mask
    if length(params.mask)>0
        % Apply mask to metric data to return
        metric_data=metric_data(params.mask);
        if params.plot || length(params.output_file)>0
            % Index of all faces
            all_faces=[1:length(display_metric_data)];
            % Get faces that aren't masked
            unmasked_faces=setdiff(all_faces,params.mask);
            % For plotting set unmasked data to extra color - gray
            display_metric_data(unmasked_faces)=numcol;
        end
    end

    % Apply threshold
    if abs(params.threshold)>0
        % Negative threshold
        if params.threshold<0
            if params.plot || length(params.output_file)>0
                % For plotting set data over threshold to extra color - gray
                display_metric_data(find(metric_data>params.threshold))=numcol;
            end
            % Apply threshold to metric data to return
            metric_data=metric_data(find(metric_data<=params.threshold));
        % Positive threshold
        else
            if params.plot || length(params.output_file)>0
                % For plotting set data under threshold to extra color - gray
                display_metric_data(find(metric_data<params.threshold))=numcol;
            end
            % Apply threshold to metric data to return
            metric_data=metric_data(find(metric_data>=params.threshold));
        end
    end

end

if params.plot || length(params.output_file)>0
    % Set face color data
    set(hp,'FaceVertexCData', display_metric_data);

    % Save plot to file
    if length(params.output_file)>0
        if params.output_format=='eps'
            figure2eps(fig, params.output_file, 10, '-opengl');
        else
            saveas(fig, params.output_file, params.output_format);
        end
    end
end
