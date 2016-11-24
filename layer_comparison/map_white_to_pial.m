function white_pial_map=map_white_to_pial(white, pial)

%white_pial_map=zeros(size(pial.vertices,1),1);
% Loop through all pial vertices
%for i=1:size(pial.vertices,1)
    % For each pial vertex, compute the distance to every other white vertex
%    pial_white_diff=white.vertices-repmat(pial.vertices(i,:),size(white.vertices,1),1);
%    pial_white_dist=sqrt(sum((pial_white_diff.^2)'));
    % find the closest white vertex to this pial vertex
%    min_dist=find(pial_white_dist==min(pial_white_dist));
%    white_pial_map(i)=min_dist(1);
%end
white_pial_map=dsearchn(white.vertices,pial.vertices);
