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
                                  overwrite = NULL) {

  #library(purrr)

  if(is.null(overwrite)) {
    overwrite <- .scrattch.io_env$overwrite
  }

  verbosity <- .scrattch.io_env$verbosity

  if(!file.exists(tome)) {
    if(verbosity == 2) {
      print(paste0(tome," doesn't exist. Creating new file."))
    }
    rhdf5::h5createFile(tome)
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

  ls <- h5ls(tome) %>%
    dplyr::mutate(full_name = ifelse(group == "/",
                                     paste0(group, name),
                                     paste(group, name, sep = "/")))

  existing_objects <- ls %>%
    dplyr::filter(group == target | full_name == target)

  if(length(existing_objects$full_name) > 0) {
    if(overwrite) {

      if(verbosity == 2) {
        print(paste0("Removing existing ", target))
      }

      purrr::walk(existing_objects$full_name,
                  function(x) {
                    suppressWarnings(
                      rhdf5::h5delete(tome, x)
                    )
                  }
      )
    } else {
      if(verbosity == 2) {
        stop(paste0(target, " already exists. Set overwrite = TRUE to replace it."))
      } else if(verbosity == 1) {
        return(FALSE)
      }

    }
  }

  # check for group structure
  write_tome_group(tome,
                   target_path)

  if(verbosity == 2) {
    print(paste0("Writing ", target))
  }

  if(store_as == "vectors") {
    purrr::walk(names(df), function(x) {
      vec_target <- paste(target,x,sep = "/")
      rhdf5::h5write(df[[x]],
                     tome,
                     vec_target,
                     level = 0)
    })

  } else if(store_as == "data.frame") {
    rhdf5::h5write(df,
                   tome,
                   target,
                   level = 0)
  }

  if(verbosity == 1) {
    return(TRUE)
  }

}

#' Generate a new group in a tome file
#'
#' @param tome Path to the target tome file
#' @param target_path The target location within the tome file
#'
write_tome_group <- function(tome,
                             target_path) {

  #library(dplyr)

  verbosity <- .scrattch.io_env$verbosity

  ls <- rhdf5::h5ls(tome)

  target_path <- sub("/$","",target_path)

  target_groups <- collapse_along(unlist(strsplit(target_path,"/"))[-1])
  if(target_groups[1] == "") {
    target_groups <- target_groups[-1]
  }
  target_groups <- paste0("/",target_groups)

  # check for group structure
  for(target_group in target_groups) {
    ls <- rhdf5::h5ls(tome) %>%
      dplyr::mutate(full_name = ifelse(group == "/",
                                       paste0(group, name),
                                       paste(group, name, sep = "/")))

    # If the group doesn't exist, create it
    if(!target_group %in% ls$full_name) {
      if(verbosity == 2) {
        print(paste0("Creating Group ",target_group))
      }
      rhdf5::h5createGroup(tome, target_group)
    } else {
      if(verbosity == 2) {
        print(paste0(target_group," already exists. Skipping."))
      }
    }
  }
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
                              overwrite = NULL) {
  #library(purrr)

  if(is.null(overwrite)) {
    overwrite <- .scrattch.io_env$overwrite
  }

  verbosity <- .scrattch.io_env$verbosity

  if(!file.exists(tome)) {
    if(verbosity == 2) {
      print(paste0(tome," doesn't exist. Creating new file."))
    }
    rhdf5::h5createFile(tome)
  }

  ls <- rhdf5::h5ls(tome) %>%
    dplyr::mutate(full_name = ifelse(group == "/",
                                     paste0(group, name),
                                     paste(group, name, sep = "/")))

  existing_objects <- ls %>%
    dplyr::filter(group == target | full_name == target)

  if(length(existing_objects$full_name) > 0) {
    if(overwrite) {
      if(verbosity == 2) {
        print(paste0("Removing existing ", target))
      }
      purrr::walk(existing_objects$full_name,
                  function(x) {
                    suppressWarnings(
                      rhdf5::h5delete(tome, x)
                    )
                  }
      )
    } else {

      if(verbosity == 2) {
        stop(paste0(target, " already exists. Set overwrite = TRUE to replace it."))
      } else if(verbosity == 1) {
        return(FALSE)
      }
    }
  }

  vec <- unlist(vec)

  if(verbosity == 2) {
    print(paste0("Writing ", target))
  }

  target_path <- sub("/$","",target)
  target_path <- sub("(/.+/).+","\\1",target_path)

  write_tome_group(tome,
                   target_path)

  rhdf5::h5write(vec,
                 tome,
                 target)

  if(verbosity == 1) {
    return(TRUE)
  }
}


#' Generalized write for individual serialized objects to a tome file
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
                                  overwrite = NULL) {

  if(is.null(overwrite)) {
    overwrite <- .scrattch.io_env$overwrite
  }

  verbosity <- .scrattch.io_env$verbosity

  if(!file.exists(tome)) {
    if(verbosity == 2) {
      print(paste0(tome," doesn't exist. Creating new file."))
    }
    rhdf5::h5createFile(tome)
  }

  ls <- rhdf5::h5ls(tome) %>%
    dplyr::mutate(full_name = ifelse(group == "/",
                                     paste0(group, name),
                                     paste(group, name, sep = "/")))

  existing_objects <- ls %>%
    dplyr::filter(group == target | full_name == target)

  if(length(existing_objects$full_name) > 0) {
    if(overwrite) {
      if(verbosity == 2) {
        print(paste0("Removing existing ", target))

      }

      purrr::walk(existing_objects$full_name,
                  function(x) {
                    suppressWarnings(
                      rhdf5::h5delete(tome, x)
                    )
                  }
      )
    } else {

      if(verbosity == 2) {
        stop(paste0(target, " already exists. Set overwrite = TRUE to replace it."))

      } else if(verbosity == 1) {
        return(FALSE)
      }


    }
  }

  ser_obj <- rawToChar(serialize(obj, NULL, ascii = TRUE))

  if(verbosity == 2) {
    print(paste0("Writing ", target))
  }

  target_path <- sub("/$","",target)
  target_path <- sub("(/.+/).+","\\1",target_path)

  write_tome_group(tome,
                   target_path)

  rhdf5::h5write(ser_obj,
                 tome,
                 target)

  if(verbosity == 1) {
    return(TRUE)
  }

}

#' Generalized write for dgCMatrix objects to a tome file
#'
#'
#' @param mat The dgCMatrix object to store
#' @param tome Path to the target tome file
#' @param target The target location within the tome file
#' @param overwrite logical, whether or not to overwrite existing objects. Default = FALSE.
#' @param compression_level The compression level to use for long vectors x and i between 0 and 9. Default = 4.
#'
write_tome_dgCMatrix <- function(mat,
                                 tome,
                                 target,
                                 overwrite = NULL,
                                 compression_level = 4) {

  #library(Matrix)

  if(is.null(overwrite)) {
    overwrite <- .scrattch.io_env$overwrite
  }

  verbosity <- .scrattch.io_env$verbosity

  target_path <- sub("/$","",target)
  write_tome_group(tome,
                   target_path)

  if(verbosity == 2) {
    print(paste0("Writing ",target,"/x."))
  }
  rhdf5::h5createDataset(tome,
                         dataset = paste0(target,"/x"),
                         dims = length(mat@x),
                         chunk = 1000,
                         level = compression_level)
  rhdf5::h5write(mat@x,
                 tome,
                 paste0(target,"/x"))

  # data indices
  if(verbosity == 2) {
    print(paste0("Writing ",target,"/i."))
  }
  rhdf5::h5createDataset(tome,
                         dataset = paste0(target,"/i"),
                         dims = length(mat@x),
                         chunk = 1000,
                         level = compression_level)
  rhdf5::h5write(mat@i,
                 tome,
                 paste0(target,"/i"))

  # data index pointers
  if(verbosity == 2) {
    print(paste0("Writing ",target,"/p."))
  }
  rhdf5::h5write(mat@p,
                 tome,
                 paste0(target,"/p"))

  rhdf5::h5write(c(nrow(mat), ncol(mat)),
                 tome,
                 paste0(target,"/dims"))
}
