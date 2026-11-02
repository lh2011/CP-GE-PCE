% n-times repeated knot insertion at intermediate knot values
function nrb = nrbrefine1D(nrb, times)
    arguments
        nrb 
        times = 1;
    end
    for itimes = 1:times
        knots = unique(nrb.knots);
        newknots = knots(2:end) - diff(knots)/2;
        nrb = nrbkntins(nrb, newknots);
    end
end