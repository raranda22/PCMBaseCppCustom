# Early-burst plus Ornstein–Uhlenbeck hybrid model (EHD)

#' @export
PCMParentClasses.EHD <- function(model) {
  c("GaussianPCM", "PCM")
}

#' @export
PCMDescribe.EHD <- function(model, ...) {
  "Early-burst plus Ornstein-Uhlenbeck hybrid model"
}

#' @export
PCMDescribeParameters.EHD <- function(model, ...) {
  list(
    X0 = "trait values at the root",
    rho = "Early burst rate parameter",
    H = "adaptation rate matrix (OU component)",
    Theta = "long-term optimum (OU component)",
    Sigma_x = "Upper triangular factor of the unit-time variance rate",
    Sigmae_x = "Upper triangular factor of the non-heritable variance or measurement error")
}

#' @export
PCMListParameterizations.EHD <- function(model, ...) {
  list(
    X0 = list(
      c("VectorParameter", "_Global"),
      c("VectorParameter", "_Fixed", "_Global"),
      c("VectorParameter", "_AllEqual", "_Global"),
      c("VectorParameter", "_Omitted")),
    rho = list(
      c("ScalarParameter", "_NonNegative"),
      c("ScalarParameter", "_NonNegative", "_Global")),
    H = list(
      c("MatrixParameter"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Global"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global")),
    Theta = list(
      c("VectorParameter"),
      c("VectorParameter", "_Global")),
    Sigma_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global")),
    Sigmae_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Omitted"))
  )
}

#' @export
PCMSpecify.EHD <- function(model, ...) {
  spec <- list(
    X0 = structure(0.0, class = c('VectorParameter', '_Global'),
                   description = 'trait values at the root'),
    rho = structure(0.0, class = c('ScalarParameter', '_NonNegative'),
                    description = 'Early burst rate parameter'),
    H = structure(0.0, class = c('MatrixParameter', '_WithNonNegativeDiagonal'),
                  description = 'adaptation rate matrix'),
    Theta = structure(0.0, class = c('VectorParameter'),
                      description = 'long-term optimum'),
    Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                        description = 'Upper triangular factor of the unit-time variance rate'),
    Sigmae_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                         description = 'Upper triangular factor of the non-heritable variance'))
  attributes(spec) <- attributes(model)
  if(is.null(names(spec))) names(spec) <- c('X0', 'rho', 'H', 'Theta', 'Sigma_x', 'Sigmae_x')
  spec
}

#' @export
PCMCond.EHD <- function(
    tree, model, r = 1, metaI = PCMInfo(NULL, tree, model, verbose = verbose),
    verbose=FALSE) {
  
  Sigma_x <- GetSigma_x(model, "Sigma", r)
  Sigma <- Sigma_x %*% t(Sigma_x)
  
  if(!is.null(model$Sigmae_x)) {
    Sigmae_x <- GetSigma_x(model, "Sigmae", r)
    Sigmae <- Sigmae_x %*% t(Sigmae_x)
  } else {
    Sigmae <- NULL
  }
  
  if(!is.null(model$rho)) {
    rho <- if(is.Global(model$rho)) model$rho else model$rho[r]
  } else {
    rho <- 0
  }
  
  if(!is.null(model$H)) {
    H <- if(is.Global(model$H)) model$H else model$H[,,r]
  } else {
    H <- matrix(0, nrow(Sigma), ncol(Sigma))
  }
  
  if(!is.null(model$Theta)) {
    Theta <- if(is.Global(model$Theta)) model$Theta else model$Theta[,r]
  } else {
    Theta <- rep(0, nrow(Sigma))
  }
  
  nHs <- nodeHeights(tree)
  
  PCMCondVOU_EHD(H = H, Theta = Theta, Sigma = Sigma, Sigmae = Sigmae, 
                 rho = rho, node_heights = nHs,
                 threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8))
}

#' @export
PCMCondVOU_EHD <- function(
    H, Theta, Sigma, Sigmae = NULL, Sigmaj = NULL, xi = NULL,
    rho = 0, node_heights = NULL,
    e_Ht = NULL,
    threshold.Lambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8)) {
  
  k <- nrow(Sigma)
  
  function(t, edgeIndex, metaI, regime) {
    if(is.null(e_Ht)) {
      if(max(abs(H)) < threshold.Lambda_ij) {
        # BM case with Early Burst scaling
        Phi <- diag(k)
        node_height <- node_heights[regime]
        eb_factor <- if(rho == 0) 1.0 else (1 - exp(-rho * node_height)) / rho
        V <- Sigma * eb_factor * t
        omega <- Theta * t
      } else {
        # OU case with Early Burst scaling
        Phi <- expm::expm(-H * t)
        node_height <- node_heights[regime] 
        eb_factor <- if(rho == 0) 1.0 else (1 - exp(-rho * node_height)) / rho
        
        # Variance integration for OU with EB scaling
        V_integral <- solve(H + t(H)) %*% (diag(k) - expm::expm(-(H + t(H)) * t))
        V <- Sigma * eb_factor * V_integral
        omega <- (diag(k) - Phi) %*% Theta
      }
    } else {
      Phi <- e_Ht[[edgeIndex]]
      node_height <- node_heights[regime]
      eb_factor <- if(rho == 0) 1.0 else (1 - exp(-rho * node_height)) / rho
      V <- Sigma * eb_factor
      omega <- (diag(k) - Phi) %*% Theta
    }
    
    if(!is.null(Sigmae)) {
      Sigmae_regime <- Sigmae
    } else {
      Sigmae_regime <- NULL
    }
    
    list(omega = omega, Phi = Phi, V = V, Sigmae = Sigmae_regime)
  }
}

#' @export
PCMInfo.EHD <- function(
    X, tree, model,
    SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
    verbose = FALSE, preorder = NULL, ...) {
  
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }
  
  res <- NextMethod()
  res$nodeHeights <- nodeHeights(tree)
  res
}
