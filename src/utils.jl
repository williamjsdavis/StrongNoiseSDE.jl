# Utility functions

# Histogram integration
# TODO: Replace this with something more efficient
function integrate(x, y)
    retval = (x[2] - x[1]) * (y[1] + y[2])
    for i in 2:(length(y) - 1)
        retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return 0.5 * retval
end

