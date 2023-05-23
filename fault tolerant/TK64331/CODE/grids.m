function gridtype = grids(r,n,xlocs,ylocs,initial,yfinal)
for i1 = 1:(n+1)
    if (xlocs(i1) >= initial) && (xlocs(i1) <= (2*r))
        if (ylocs(i1) >= initial) && (ylocs(i1) <= (2*r))
            gridtype(i1) = 1;
        elseif (ylocs(i1) > (2*r)) && (ylocs(i1) <= (4*r))
            gridtype(i1) = 2;
        elseif (ylocs(i1) > (4*r)) && (ylocs(i1) <= yfinal)
            gridtype(i1) = 3;
        end
    elseif (xlocs(i1) > (2*r)) && (xlocs(i1) <= (4*r))
        if (ylocs(i1) > initial) && (ylocs(i1) <= (2*r))
            gridtype(i1) = 4;
        elseif (ylocs(i1) > (2*r)) && (ylocs(i1) <= (4*r))
            gridtype(i1) = 5;
        elseif (ylocs(i1) > (4*r)) && (ylocs(i1) <= yfinal)
            gridtype(i1) = 6;
        end
    elseif (xlocs(i1) > (4*r)) && (xlocs(i1) <= (6*r))
        if (ylocs(i1) >= initial) && (ylocs(i1) <= (2*r))
            gridtype(i1) = 7;
        elseif (ylocs(i1) > (2*r)) && (ylocs(i1) <= (4*r))
            gridtype(i1) = 8;
        elseif (ylocs(i1) > (4*r)) && (ylocs(i1) <= yfinal)
            gridtype(i1) = 9;
        end
    elseif (xlocs(i1) > (6*r)) && (xlocs(i1) <= (8*r))
        if (ylocs(i1) >= initial) && (ylocs(i1) <= (2*r))
            gridtype(i1) = 10;
        elseif (ylocs(i1) > (2*r)) && (ylocs(i1) <= (4*r))
            gridtype(i1) = 11;
        elseif (ylocs(i1) > (4*r)) && (ylocs(i1) <= yfinal)
            gridtype(i1) = 12;
        end
    elseif (xlocs(i1) > (8*r)) && (xlocs(i1) <= (10*r))
        if (ylocs(i1) >= initial) && (ylocs(i1) <= (2*r))
            gridtype(i1) = 13;
        elseif (ylocs(i1) > (2*r)) && (ylocs(i1) <= (4*r))
            gridtype(i1) = 14;
        elseif (ylocs(i1) > (4*r)) && (ylocs(i1) <= yfinal)
            gridtype(i1) = 15;
        end
    else
        if (ylocs(i1) > initial) && (ylocs(i1) <= (2*r))
            gridtype(i1) = 16;
        elseif (ylocs(i1) > (2*r)) && (ylocs(i1) <= (4*r))
            gridtype(i1) = 17;
        elseif (ylocs(i1) > (4*r)) && (ylocs(i1) <= yfinal)
            gridtype(i1) = 18;
        end
    end
end
end