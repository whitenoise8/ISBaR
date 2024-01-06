#' file: utils_blm_reg.R
#' author: Cristian Castiglione
#' creation: 05/01/2023
#' last change: 05/01/2023

# Function for initializing the prior parameters of the BLM
init.blm.prior = function (prior) {
  default = list(
    a = 2.01,
    b = 1.01,
    v = 0.1)
  
  if (class(prior) != "list") {
    # Print a warning if `control` is not a list
    warning("The `prior` parameter must be a list \n",
            "Default parameters will be used for the execution.")
  } else {
    if (length(prior) != 0) {
      # Get the default and prior attributes
      nms.default = names(default)
      nms.prior = names(prior)
      nms.undef = nms.prior[!nms.prior %in% nms.default]
      
      # Set the custom hyperparameters
      default[names(prior)] = prior
      
      # Print a warning if some of the prior parameters are not allowed 
      if (length(nms.undef)) {
        warning("Unknown names in prior parameters: ", 
                paste(nms.undef, collapse = ", "))
      }
    }
  }
  
  # Return the control parameters
  return (default)
}

# Function for initializing the control parameters of the algorithms
init.blm.control = function (control) {
  
  # Default control parameters
  default = list(
    niter = 5000,
    burn = 2500,
    thin = 1,
    verbose = TRUE,
    report = 500)
  
  if (class(control) != "list") {
    # Print a warning if `control` is not a list
    warning("The `control` parameter must be a list \n",
            "Default parameters will be used for the execution.")
  } else {
    if (length(control) != 0) {
      # Get the default and control attributes
      nms.default = names(default)
      nms.control = names(control)
      nms.undef = nms.control[!nms.control %in% nms.default]
      
      # Set the custom parameters
      default[names(control)] = control
      
      # Print a warning if some of the control parameters are not allowed 
      if (length(nms.undef)) {
        warning("Unknown names in control: ", 
                paste(nms.undef, collapse = ", "))
      }
    }
  }
  
  # Return the control parameters
  return (default)
}
