function [ index_centres ] = select_init_centers( nbElements, nbIndices )
    rearrange_positions = randperm(nbElements);
    index_centres = rearrange_positions(1:nbIndices);
end

