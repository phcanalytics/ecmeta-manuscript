my_update <- function (object, call, ...) {
  extras <- match.call(expand.dots = FALSE)$...
  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    ## do these individually to allow NULL to remove entries.
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, parent.frame())
}

#' Update treatment effect estimates
#' 
#' Update models used to estimate average treatment effects (ATEs) with survival data.
#' This will re-estimate the function used to create `object` and is helpful for
#' performing sensitivity analyses where only a few arguments need to be modified. It
#' does this by extracting the original function "call". 
#' 
#' @param object An object of the appropriate class. 
#' @return An updated version of `object`.
#' @name update_ate
#' @export
update.surv_ate <- function (object, ...) {
  
  call <- attr(object, "call")
  my_update(object, call, ...)
}

#' @rdname update_ate
#' @export
update.grouped_surv_ate <- function (object, ...) {
  update.surv_ate(object, ...)
}

#' @rdname update_ate
#' @export
update.grouped_ps_surv <- function (object, ...) {
  call <- object$call
  my_update(object, call, ...)
}