# created on Dec. 15, 2015
#  (1) obtain the estimate of corrected AUC and its confidence interval
#       based on Reiser's (2000) method (balanced case with 2 replicates per
#       subject)
#
# datFrame - data frame with at least 2 columns
#   'y' -- observations
#   'subjID' -- subject id
#   'grp' -- group id: 1 - case; 0 - control;
#   'myrep' -- replication id: 1, 2, 3...

AUCest.Reiser <- function (datFrame, 
    sidVar = "subjID", 
    obsVar = "y", 
    grpVar = "grp", 
    repVar = "myrep", alpha = 0.05) 
{
    idVec = datFrame[, c(sidVar)]
    myrep = datFrame[, c(repVar)]
    datFrame = datFrame[order(idVec, myrep), ]
    idVec = datFrame[, c(sidVar)]
    myrep = datFrame[, c(repVar)]
    u.rep = sort(unique(myrep))
    n.rep = length(u.rep)
    if (prod(c(1:n.rep) == u.rep) != TRUE) {
        stop("replication indicator should be consecutive integer starting from 1\n!")
    }
    mygrp = datFrame[, c(grpVar)]
    u.grp = sort(unique(mygrp))
    n.grp = length(u.grp)
    if (prod(c(0, 1) == u.grp) != TRUE) {
        stop("group indicator should be 0 or 1\n!")
    }
    myscore <- datFrame[, c(obsVar)]
    scoreMat = myscore[myrep == u.rep[1]]
    for (i in 2:n.rep) {
        scoreMat = cbind(scoreMat, myscore[myrep == u.rep[i]])
    }
    m.score = apply(scoreMat, 1, mean, na.rm = TRUE)
    grp1 = datFrame[which(myrep == 1), c(grpVar)]
    Xidotbar = m.score[which(grp1 == 1)]
    Yidotbar = m.score[which(grp1 == 0)]
    X <- myscore[mygrp == 1]
    Y <- myscore[mygrp == 0]
    muX <- mean(X, na.rm = TRUE)
    muY <- mean(Y, na.rm = TRUE)
    nX = sum(grp1 == 1, na.rm = TRUE)
    nY = sum(grp1 == 0, na.rm = TRUE)
    omega2 = (sum((Xidotbar - muX)^2, na.rm = TRUE) + sum((Yidotbar - 
        muY)^2, na.rm = TRUE))/(nX + nY - 2)
    scoreMatX = scoreMat[which(grp1 == 1), ]
    scoreMatY = scoreMat[which(grp1 == 0), ]
    p = n.rep
    tau2 = (sum((scoreMatX[, 1] - Xidotbar)^2, na.rm = TRUE) + 
        sum((scoreMatX[, 2] - Xidotbar)^2, na.rm = TRUE) + sum((scoreMatY[, 
        1] - Yidotbar)^2, na.rm = TRUE) + sum((scoreMatY[, 2] - 
        Yidotbar)^2, na.rm = TRUE))/((nX + nY) * (p - 1))
    sigma2 = omega2 - tau2/p
    if (sigma2 > 0) {
        mu = muX - muY
        part1 = omega2 * (1/nX + 1/nY)/(2 * sigma2)
        part2.1 = 2 * omega2^2/(nX + nY - 2)
        part2.2 = 2 * tau2^2/(p^2 * (nX + nY) * (p - 1))
        sd.delta = sqrt(part1 + mu^2 * (part2.1 + part2.2)/(8 * 
            sigma2^3))
        delta = mu/sqrt(sigma2 * 2)
        za = qnorm(1 - alpha/2)
        CI.low.delta = delta - za * sd.delta
        CI.upp.delta = delta + za * sd.delta
        AUC.hat = pnorm(delta)
        CI.low = pnorm(CI.low.delta)
        CI.upp = pnorm(CI.upp.delta)
        sd.AUC.hat = dnorm(delta) * sd.delta
        res <- list(AUC.c = AUC.hat, sd.AUC.c = sd.AUC.hat, AUC.c.low = CI.low, 
            AUC.c.upp = CI.upp)
    }
    else {
        res <- list(AUC.c = NA, sd.AUC.c = NA, AUC.c.low = NA, AUC.c.upp = NA)
    }
    invisible(res)
}
