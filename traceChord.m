% hmm forgot i wrote this. i guess its redundant with extendRegularHexChord
% now. but that one is more congruent with the other extend chord things so
% ill keep that and deprecate this.
function chord = traceChord(data, fseed, hexseed, stopOnSingular)
    chord = false(size(data.H,1),1);
    
    currhex = hexseed;
    currface = fseed;
    loopflag = false;
    while true
        chord(currhex) = true;
        
        nextface = data.H2Farray(currhex,data.H2F_flip(data.H2F6(currhex, currface)));
        if data.isBoundaryFace(nextface)
            break;
        end
        nexthex = setdiff(data.F2Harray(nextface,[1 3]), currhex);
        
        if chord(nexthex) & loopflag
            % already in chord and paralell.
            break;
        elseif chord(nexthex) 
            % already in chord. could be perpendicular?
            loopflag = true;
        else
            loopflag = false;
        end
        
        if stopOnSingular
            if any(data.isSingularVertex(data.H(nexthex,:)))
                chord(nexthex) = true;
                break;
            end
        end
        
        currhex = nexthex;
        currface = nextface;
    end
    
    
end