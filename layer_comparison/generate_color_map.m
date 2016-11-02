function color_vector=generate_color_map(numColors, control_pts)

    color_vector = [];

    
    for color=1:numColors
        cx = 10 + ((color-1)/(numColors-1))*(440);
        for j=1:size(control_pts,2)-1
            thisPt = control_pts(:,j);
            nextPt = control_pts(:,j+1);
            thisPtX=10 + (460-20) / (size(control_pts,2) - 1) * (j-1);
            nextPtX=10 + (460-20) / (size(control_pts,2) - 1) * (j);
                
            % Find the points to interpolate between
            if(cx >= thisPtX && cx <= nextPtX)
                cyInterp = thisPt + (cx - thisPtX)/(nextPtX - thisPtX)*(nextPt - thisPt);

                if(cyInterp > 350 - 10)
                    cyInterp = 350 - 10;
                elseif (cyInterp < 10)
                    cyInterp = 10;
                end
                color_vector(end+1,:) = round(255*(1-(cyInterp-10)/(350-20)));
                break

            else
                continue            
            end
        end
    end

