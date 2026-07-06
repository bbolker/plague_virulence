##' Update a metapop_run by modifying arguments and re-running
##'
##' Retrieves the \code{call} attribute stored on a \code{metapop_run} object,
##' substitutes any named arguments supplied via \code{...}, and re-evaluates
##' the modified call in the caller's environment.  Unrecognised argument names
##' (not present in the original call) are appended rather than replacing an
##' existing entry.
##'
##' @param object a \code{metapop_run} object (output of \code{discrete_run()})
##' @param ... named arguments to replace (or add) in the stored call
##' @param evaluate logical; if \code{TRUE} (default) the updated call is
##'   evaluated and its result returned; if \code{FALSE} the modified call
##'   object is returned without evaluating
##' @return a new \code{metapop_run} object, or (if \code{evaluate = FALSE})
##'   the updated \code{call} object
##' @export
update.metapop_run <- function(object, ..., evaluate = TRUE) {
  cc <- attr(object, "call")
  if (is.null(cc))
    stop("'object' has no 'call' attribute")
  extras <- match.call(expand.dots = FALSE)$`...`
  if (length(extras) > 0L) {
    existing <- names(extras) %in% names(cc)
    for (nm in names(extras)[existing])
      cc[[nm]] <- extras[[nm]]
    if (any(!existing))
      cc <- as.call(c(as.list(cc), extras[!existing]))
  }
  if (evaluate) eval(cc, parent.frame()) else cc
}
