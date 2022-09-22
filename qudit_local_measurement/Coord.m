function [Coord_psi] = Coord(psi)
    k  = 0;
    dimension = length(psi);
    for j = 1 : dimension
        if psi(j) ~= 0
            k = k + 1;
            Coord_psi(k,1) = j-1; %%s: contains all of elements with non-zero values, start from 0-->0...0
            Coord_psi(k,2) = psi(j);
        end
    end
end
