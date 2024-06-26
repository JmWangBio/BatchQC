#' This function allows you to add transformed count matrix to the SE object
#' @param se SummarizedExperiment Object
#' @param method Transformation Method, currently logit only
#' @param assay_to_transform Which SE assay to do transformation on
#' @param output_assay_name name for the resulting transformed assay
#' @return the original SE object with transformed assay appended
#' @import SummarizedExperiment
#'
#' @export
transform_SE <- function(se, method, assay_to_transform, output_assay_name) {
    se <- se
    if (method == 'logit') {
        assays(se)[[output_assay_name]] <-
            log(assays(se)[[assay_to_transform]] / (1 - assays(se)[[assay_to_transform]]))
    }
    return(se)
}
