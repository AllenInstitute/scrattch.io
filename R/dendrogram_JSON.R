#' Convert a dendrogram object to a nested list
#'
#' This function preserves node and leaf attributes as data.frame objects by converting each
#' nested dendrogram into a list containing the attribute data.frame and a nested list of children.
#'
#' This provides a structure suitable for export to JSON.
#'
#' @param dend a dendrogram object
#'
#' @return a list object, as described above.
#' @export
dend_to_list <- function(dend) {

  node_attributes <- as.data.frame(attributes(dend)[names(attributes(dend)) != "class"])
  node_attributes <- unique(node_attributes[,names(node_attributes) != "names"])
  #print(node_attributes)
  if("leaf" %in% names(node_attributes)) {
    #print("leaf")
    return(list(leaf_attributes = node_attributes))
  } else {
    #print("node")
    y <- dend
    attributes(y) <- NULL
    class(y) <- "list"
    children <- y

    dend <- list(node_attributes = node_attributes,
                 children = children)

    if(length(dend$children) > 1) {
      for(i in 1:length(dend$children)) {
        dend$children[[i]] <- dend_to_list(dend$children[[i]])
      }
    }
    return(dend)
  }

}


#' Convert a dendrogram into a nested JSON object
#'
#' @param dend a dendrogram object
#'
#' @return a JSON text object
#' @export
#'
dend_to_JSON <- function(dend) {
  dend_list <- dend_to_list(dend)
  jsonlite::toJSON(dend_list,
                   complex = "list",
                   pretty = TRUE)
}
