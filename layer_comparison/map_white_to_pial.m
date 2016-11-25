function white_pial_map=map_white_to_pial(white, pial)

white_pial_map=dsearchn(white.vertices,pial.vertices);
