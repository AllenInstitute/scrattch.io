#' Generalized write for data.frames to a tome file
#'
#' This function currently only works in an overwrite mode. Anything at the target
#' location will be removed, and the new df will be written.
#'
#' @param df The data.frame to store
#' @param tome Path to the target tome file
#' @param target The target location within the tome file
#' @param store_as Either "data.frame", which will store as a compound object; or
#' "vectors", which stores each column as a separate object.
#' @param overwrite logical, whether or not to overwrite existing objects. Default = FALSE.
#'
write_tome_data.frame <- function(df,
                                  tome,
                                  target,
                                  store_as = "vectors",
                                  overwrite = FALSE) {

  library(rhdf5)
  library(purrr)
  library(h5)

  H5close()

  if(!file.exists(tome)) {
    print(paste0(tome," doesn't exist. Creating new file."))
    h5createFile(tome)
    H5close()
  }

  target_path <- sub("/$","",target)

  if(store_as == "vectors") {
    target_path <- target_path
  } else if(store_as == "data.frame") {
    target_path <- sub("(/.+/).+","\\1",target_path)
  }

  target_path <- sub("/$","",target_path)

  target_groups <- collapse_along(unlist(strsplit(target_path,"/"))[-1])
  target_groups <- paste0("/",target_groups)

  ## Now, need to check for existing vectors (for vector write)
  ## or existing frame (for data.frame write)
  #
  #   existing_anno <- ls %>% filter(group == "/sample_meta/anno")
  #   anno_cols <- names(anno)

  # If overwrite, remove everything
  #if(mode == "overwrite") {

  ls <- h5ls(tome) %>%
    mutate(full_name = ifelse(group == "/",
                              paste0(group, name),
                              paste(group, name, sep = "/")))

  existing_objects <- ls %>%
    filter(group == target | full_name == target)

  if(length(existing_objects$full_name) > 0) {
    if(overwrite) {

      print(paste0("Removing existing ", target))

      walk(existing_objects$full_name,
           function(x) {
             suppressWarnings(
               h5_delete(tome, x)
             )
           }
      )
    } else {

      stop(paste0(target, " already exists. Set overwrite = TRUE to replace it."))

    }
  }

  # check for group structure
  for(target_group in target_groups) {
    ls <- h5ls(tome) %>%
      mutate(full_name = ifelse(group == "/",
                                paste0(group, name),
                                paste(group, name, sep = "/")))

    if(!target_group %in% ls$full_name) {
      print(paste0("Creating Group ",target_group))
      h5createGroup(tome, target_group)
    }
  }

  print(paste0("Writing ", target))

  if(store_as == "vectors") {
    walk(names(df), function(x) {
      vec_target <- paste(target,x,sep = "/")
      h5write(df[[x]],
              tome,
              vec_target)
    })

  } else if(store_as == "data.frame") {
    h5write(df,
            tome,
            target)
  }

  # }
  # alternate modes for future use.
  # else if(mode == "add") {
  #
  #   if(nrow(existing_anno) > 0) {
  #     anno_cols <- anno_cols[!anno_cols %in% existing_anno$name]
  #
  #     if(length(anno_cols) > 0) {
  #       print(paste0("Adding columns ", paste(anno_cols, collapse = ", "), " to /sample_meta/anno."))
  #
  #       # Order to match existing annotations
  #       existing_sample_names <- h5read(tome, "/sample_meta/anno/sample_name")
  #
  #       anno <- anno %>%
  #         select(one_of("sample_name", anno_cols))
  #
  #       anno <- anno[match(anno$sample_name, existing_sample_names),]
  #
  #       # Write the new columns
  #       walk(anno_cols, function(x) {
  #         target <- paste0("/sample_meta/anno/",x)
  #         h5write(anno[[x]],
  #                 tome,
  #                 target)
  #       })
  #     } else {
  #       print("No new columns to add. If you want to replace values, try mode = 'replace' to replace the whole anno table.")
  #     }
  #   } else {
  #     print("Writing /sample_meta/anno")
  #     anno_cols <- names(anno)
  #     walk(anno_cols, function(x) {
  #       target <- paste0("/sample_meta/anno/",x)
  #       h5write(anno[[x]],
  #               tome,
  #               target)
  #     })
  #   }
  # }

}

#' Generalized write for individual vector objects to a tome file
#'
#' @param vec The vector to store
#' @param tome Path to the target tome file
#' @param target The target location within the tome file
#' @param overwrite logical, whether or not to overwrite existing objects. Default = FALSE.
#'
write_tome_vector <- function(vec,
                              tome,
                              target,
                              overwrite = FALSE) {
  library(rhdf5)
  library(purrr)
  library(h5)

  H5close()

  if(!file.exists(tome)) {
    print(paste0(tome," doesn't exist. Creating new file."))
    h5createFile(tome)
    H5close()
  }

  ls <- h5ls(tome) %>%
    mutate(full_name = ifelse(group == "/",
                              paste0(group, name),
                              paste(group, name, sep = "/")))

  existing_objects <- ls %>%
    filter(group == target | full_name == target)

  if(length(existing_objects$full_name) > 0) {
    if(overwrite) {

      print(paste0("Removing existing ", target))

      walk(existing_objects$full_name,
           function(x) {
             suppressWarnings(
               h5_delete(tome, x)
             )
           }
      )
    } else {

      stop(paste0(target, " already exists. Set overwrite = TRUE to replace it."))

    }
  }

  vec <- unlist(vec)

  print(paste0("Writing ", target))

  h5write(vec,
          tome,
          target)

}


#' Generalized write for individual vector objects to a tome file
#'
#' Useful for R objects that can't easily be coerced to a data.frame or a vector, like lists or S3 classes.
#' Think of it like saveRDS() for HDF5 files.
#'
#' @param obj The R object to store
#' @param tome Path to the target tome file
#' @param target The target location within the tome file
#' @param overwrite logical, whether or not to overwrite existing objects. Default = FALSE.
#'
write_tome_serialized <- function(obj,
                                  tome,
                                  target,
                                  overwrite = FALSE) {

  library(rhdf5)
  library(purrr)
  library(h5)

  H5close()

  if(!file.exists(tome)) {
    print(paste0(tome," doesn't exist. Creating new file."))
    h5createFile(tome)
    H5close()
  }

  ls <- h5ls(tome) %>%
    mutate(full_name = ifelse(group == "/",
                              paste0(group, name),
                              paste(group, name, sep = "/")))

  existing_objects <- ls %>%
    filter(group == target | full_name == target)

  if(length(existing_objects$full_name) > 0) {
    if(overwrite) {

      print(paste0("Removing existing ", target))

      walk(existing_objects$full_name,
           function(x) {
             suppressWarnings(
               h5_delete(tome, x)
             )
           }
      )
    } else {

      stop(paste0(target, " already exists. Set overwrite = TRUE to replace it."))

    }
  }

  ser_obj <- rawToChar(serialize(obj, ascii = TRUE))

  print(paste0("Writing ", target))

  h5write(ser_obj,
          tome,
          target)

}
