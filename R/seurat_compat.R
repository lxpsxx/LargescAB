get_scab_assay_layer <- function(Obejct, assay = "RNA", layer = "data") {
  assay_obj <- Obejct[[assay]]

  if (is.null(assay_obj)) {
    stop(sprintf("Assay '%s' is not found in the input Seurat object.", assay))
  }

  if (inherits(assay_obj, "Assay5")) {
    return(SeuratObject::LayerData(assay_obj, layer = layer))
  }

  if (!(layer %in% methods::slotNames(assay_obj))) {
    stop(sprintf("Layer/slot '%s' is not available in assay '%s'.", layer, assay))
  }

  methods::slot(assay_obj, layer)
}


get_scab_graph <- function(Obejct, graph_name = "RNA_snn") {
  if (!(graph_name %in% names(Obejct@graphs)) || is.null(Obejct@graphs[[graph_name]])) {
    stop(sprintf("Graph '%s' is not found, please run FindNeighbors first.", graph_name))
  }

  as.matrix(Obejct@graphs[[graph_name]])
}


get_scab_graph_sparse <- function(Obejct, graph_name = "RNA_snn") {
  if (!(graph_name %in% names(Obejct@graphs)) || is.null(Obejct@graphs[[graph_name]])) {
    stop(sprintf("Graph '%s' is not found, please run FindNeighbors first.", graph_name))
  }

  methods::as(Matrix::Matrix(Obejct@graphs[[graph_name]], sparse = TRUE), "dgCMatrix")
}
