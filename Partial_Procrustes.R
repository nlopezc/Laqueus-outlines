#Partial Procrustes based on Momocs' code for GPA

.vecs_param <- function(r1, i1, r2, i2) {
  x <- c(r1, i1, r2, i2)
  if (!is.numeric(x)) {
    stop("4 numeric must be passed")
  }
  if (length(x) != 4) {
    stop("4 numeric must be passed")
  }
  r.norms <- sqrt((r2^2 + i2^2))/sqrt((r1^2 + i1^2))
  d1 <- sqrt(sum(r1^2 + i1^2))
  d2 <- sqrt(sum(r2^2 + i2^2))
  return(list(r.norms = d1/d2, d.angle = atan2(i2, r2) - atan2(i1,
                                                               r1)))
}

part_procrustes_default <- function(x, tol = 1e-05,  coo=NULL) {
  if (is.list(x))
    A <- l2a(x)
  else
    A <- x
  A <- ldk_check(A)
  # directly borrowed from Claude
  p <- dim(A)[1]
  k <- dim(A)[2]
  n <- dim(A)[3]
  if (p <= 2)
    stop("fgProcrustes makes sense with at least 3 points")
  # we prepare an array to save results
  temp2 <- temp1 <- array(NA, dim = c(p, k, n))
  Siz <- numeric(n)
  for (i in 1:n) {
    Siz[i] <- coo_centsize(A[, , i])
    temp1[, , i] <- coo_center(A[, , i])
  }
  sf <- NA
  M <- temp1[, , 1]
  # we do the procrustes alignment on every shape
  for (i in 1:n) {
    temp1[, , i] <- fProcrustes(temp1[, , i], M)$coo1
  }
  M <- mshapes(temp1)
  Qm1 <- dist(t(matrix(temp1, k * p, n)))
  Qd <- Qi <- Q <- sum(Qm1)
  # we initialize the counter
  iter <- 0
  sc <- rep(1, n)
  # and we loop
  while (abs(Q) > tol) {
    for (i in 1:n) {
      Z1 <- temp1[, , i]
      sv <- svd(t(M) %*% Z1)
      U <- sv$v
      V <- sv$u
      Delt <- sv$d
      sig <- sign(det(t(Z1) %*% M))
      Delt[k] <- sig * abs(Delt[k])
      V[, k] <- sig * V[, k]
      phi <- U %*% t(V)
      beta <- sum(Delt)
      temp1[, , i] <- X <- sc[i] * Z1 %*% phi
    }
    M <- mshapes(temp1)
    for (i in 1:n) {
      sf[i] <- sqrt(
        sum(diag(temp1[, , i] %*% t(M))) / (sum(diag(M %*% t(M))) *
                                              sum(diag(temp1[, , i] %*% t(temp1[, , i])))))
      temp2[, , i] <- sf[i] * temp1[, , i]
    }
    M <- mshapes(temp2)
    sc <- sf * sc
    Qm2 <- dist(t(matrix(temp2, k * p, n)))
    Qd[iter] <- Q <- sum(Qm1) - sum(Qm2)
    Qm1 <- Qm2
    Qi[iter] <- sum(Qm2)
    iter <- iter + 1
    message("iteration: ", iter, "\tgain:", signif(abs(Q), 5))
    temp1 <- temp2
  } # end of the big loop
  list(rotated = temp2,
       iterationnumber = iter, Q = Q, Qi = Qi,
       Qd = Qd,
       intereuclidean.dist = Qm2,
       mshape = coo_centsize(mshapes(temp2)),
       cent.size = Siz)
}

part_procrustes_out <- function(x, tol = 1e-10, coo=FALSE) {
  Coo <- validate(x)
  # if no $ldk defined, we convert Out into a Ldk and then
  # perform the fgProcrustes and return back an Out object.
  # THE ABOVE HAS BEEN BROKEN!

  # case where coo=FALSE and we work on the ldk THIS IS THE ONLY WORKING VERSION
  Coo2 <- coo_center(Coo)
  ref <- l2a(get_ldk(Coo2))
  nb_ldk <- dim(ref)[1]
  # case with one ldk
  if (nb_ldk == 1)
    stop("cannot apply fgProcrustes on less than three landmarks")
  # case with two ldks
  if (nb_ldk == 2) {
    message("cannot apply fgProcrustes on less than three landmarks. coo_bookstein is returned")
    return(coo_bookstein(Coo))
  }

  tar <- part_procrustes_default(ref, tol = tol)$rotated
  # would benefit to be handled by coo_baseline ?
  for (i in 1:length(Coo2)) {
    tari <- tar[, , i]
    refi <- ref[, , i]
    t1x <- tari[1, 1]
    t1y <- tari[1, 2]
    t2x <- tari[2, 1]
    t2y <- tari[2, 2]
    r1x <- refi[1, 1]
    r1y <- refi[1, 2]
    r2x <- refi[2, 1]
    r2y <- refi[2, 2]
    # translation
    t <- tari[1, ] - refi[1, ]
    refi <- coo_trans(refi, t[1], t[2])
    # rotation
    tx <- t2x - t1x
    ty <- t2y - t1y
    rx <- r2x - r1x
    ry <- r2y - r1y
    vi <- .vecs_param(rx, ry, tx, ty)
    coo_i <- Coo2$coo[[i]]
    coo_i <- coo_trans(coo_i, t[1] - t1x, t[2] - t1y)
    #coo_i <- coo_i/vi$r.norms # this line looks suspicious
    coo_i <- coo_rotate(coo_i, vi$d.angle)
    coo_i <- coo_trans(coo_i, t1x, t1y)
    Coo2$coo[[i]] <- coo_center(coo_i)
  }
  return(Coo2)
}
