iterateObs <- function(sub_index, type, cutoff, ...) {
    if(type == "spatial" & balanced_pnl) {
        
        sub_dt <- dt[time == sub_index]
        n1 <- nrow(sub_dt)
        if(n1 > 1000 & verbose){message(paste("Starting on sub index:", sub_index))}

        X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
        e <- sub_dt[, e]

        XeeXhs <- Bal_XeeXhC(d, X, e, n1, k)

    } else if(type == "spatial" & !balanced_pnl) {
        
        sub_dt <- dt[time == sub_index]
        n1 <- nrow(sub_dt)
        if(n1 > 1000 & verbose){message(paste("Starting on sub index:", sub_index))}

        X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
        e <- sub_dt[, e]
        lat <- sub_dt[, lat]; lon <- sub_dt[, lon]

        # If n1 >= 50k obs, then avoiding construction of distance matrix.
        # This requires more operations, but is less memory intensive.
        if(n1 < 5 * 10^4) {
            XeeXhs <- XeeXhC(cbind(lat, lon), cutoff, X, e, n1, k,
                kernel, dist_fn)
        } else {
            XeeXhs <- XeeXhC_Lg(cbind(lat, lon), cutoff, X, e, n1, k,
                kernel, dist_fn)
        }
    
    } else if(type == "serial") {
        sub_dt <- dt[unit == sub_index]
        n1 <- nrow(sub_dt)
        if(n1 > 1000 & verbose){message(paste("Starting on sub index:", sub_index))}

        X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE] )
        e <- sub_dt[, e]
        times <- sub_dt[, time]

        XeeXhs <- TimeDist(times, cutoff, X, e, n1, k)
    }

    XeeXhs
}

