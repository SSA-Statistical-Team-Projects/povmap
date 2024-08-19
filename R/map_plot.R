# External documentation -------------------------------------------------------

#' Visualizes Regional Disaggregated Estimates on a Map
#'
#' Function \code{map_plot} creates spatial visualizations of the estimates
#' obtained by small area estimation methods or direct estimation.
#'
#' @param object an object of type povmap, containing the estimates to be
#' visualized.
#' @param indicator optional character vector that selects which indicators
#' shall be returned: (i) all calculated indicators ("all");
#' (ii) each indicator name: "Mean", "Quantile_10", "Quantile_25", "Median",
#' "Quantile_75", "Quantile_90", "Head_Count", "Poverty_Gap", "Gini",
#' "Quintile_Share" or the function name/s of "custom_indicator/s";
#' (iii) groups of indicators: "Quantiles", "Poverty" or
#' "Inequality". Note, additional custom indicators can be
#' defined as argument for model-based approaches (see also \code{\link{ebp}})
#' and do not appear in groups of indicators even though these might belong to
#' one of the groups. If the \code{model} argument is of type "fh",
#' indicator can be set to "all", "Direct", FH", or "FH_Bench" (if povmap
#' object is overwritten by function benchmark). Defaults to "all".
#' @param MSE optional logical. If \code{TRUE}, the MSE is also visualized.
#' Defaults to \code{FALSE}.
#' @param CV optional logical. If \code{TRUE}, the CV is also visualized.
#' Defaults to \code{FALSE}.
#' @param map_obj an \code{sf, data.frame} object as defined by the
#' \pkg{sf} package on which the data should be visualized. The typical example
#' is polygon shapefile object
#' @param map_dom_id a character string containing the name of a variable in
#' \code{map_obj} that indicates the domains.
#' @param map_tab a \code{data.frame} object with two columns that match the
#' domain variable from the census data set (first column) with the domain
#' variable in the map_obj (second column). This should only be used if the IDs
#' in both objects (`map_obj` and `object`) differ.
#' @param color a \code{vector} of length 2 defining the lowest and highest
#' color in the plots.
#' @param viridis_option a one-character string between "A" and "H" requesting a 
#' viridis palette with the specified option. Defaults to NULL 
#' @param scale_points a structure defining the lowest and the highest
#' value of the colorscale. If a numeric vector of length two is given, this
#' scale will be used for every plot.
#' @param guide character passed to
#' \code{scale_colour_gradient} from \pkg{ggplot2}.
#' Possible values are "none", "colourbar", and "legend".
#' @param na.color a character string indicating the color used to plot regions with no data, 
#' passed to the na.value argument in scale_fill_gradient in ggplot2. Defaults to "transparent".  
#' @param return_data if set to \code{TRUE}, a fortified data frame including
#' the map data as well as the chosen indicators is returned. Customized maps
#' can easily be obtained from this data frame via the package \pkg{ggplot2}.
#' Defaults to \code{FALSE}.
#' @param save_file a character string that specified the prefix of the file name to save the plot. 
#' The suffix is the name of the indicator. Defaults to NULL, in which case plots are not saved. 
#' @param save_format a character string that is passed to ggsave's device argument to determine the 
#' format that the plot is saved in. Defaults to "pdf"
#' @return Creates the plots demanded, and, if selected, a fortified data.frame
#' containing the mapdata and chosen indicators.
#' @seealso \code{\link{direct}}, \code{\link{ebp}}, \code{\link{fh}},
#' \code{\link{povmapObject}}
#' @examples
#' \donttest{
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' # Generate povmap object with additional indicators; here via function ebp()
#' povmap_model <- ebp(
#'   fixed = eqIncome ~ gender + eqsize + cash +
#'     self_empl + unempl_ben + age_ben + surv_ben + sick_ben +
#'     dis_ben + rent + fam_allow + house_allow + cap_inv +
#'     tax_adj, pop_data = eusilcA_pop,
#'   pop_domains = "district", smp_data = eusilcA_smp,
#'   smp_domains = "district", threshold = 11064.82,
#'   transformation = "box.cox", L = 50, MSE = TRUE, B = 50
#' )
#'
#' # Load shape file
#' load_shapeaustria()
#'
#' # Create map plot for mean indicator - point and MSE estimates but no CV
#' map_plot(
#'   object = povmap_model, MSE = TRUE, CV = FALSE,
#'   map_obj = shape_austria_dis, indicator = c("Mean"),
#'   map_dom_id = "PB"
#' )
#'
#' # Create a suitable mapping table to use numerical identifiers of the shape
#' # file
#'
#' # First find the right order
#' dom_ord <- match(shape_austria_dis$PB, povmap_model$ind$Domain)
#'
#' # Create the mapping table based on the order obtained above
#' map_tab <- data.frame(
#'   pop_data_id = povmap_model$ind$Domain[dom_ord],
#'   shape_id = shape_austria_dis$BKZ
#' )
#'
#' # Create map plot for mean indicator - point and CV estimates but no MSE
#' # using the numerical domain identifiers of the shape file
#'
#' map_plot(
#'   object = povmap_model, MSE = FALSE, CV = TRUE,
#'   map_obj = shape_austria_dis, indicator = c("Mean"),
#'   map_dom_id = "BKZ", map_tab = map_tab
#' )
#' }
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes geom_polygon facet_wrap fortify coord_equal labs
#' @importFrom ggplot2 theme element_blank guides scale_fill_gradient
#' @importFrom ggplot2 scale_colour_gradient geom_sf scale_colour_manual guide_legend guide_colourbar
#' @importFrom viridis scale_fill_viridis
#' 

map_plot <- function(object,
                     indicator = "all",
                     MSE = FALSE,
                     CV = FALSE,
                     map_obj = NULL,
                     map_dom_id = NULL,
                     map_tab = NULL,
                     color = c("white", "red4"),
                     viridis_option = NULL, 
                     scale_points = NULL,
                     na.color = "transparent",
                     guide = "colourbar",
                     return_data = FALSE,
                     save_file = NULL,
                     save_format = "pdf") {
  if (is.null(map_obj)) {
    message(strwrap(prefix = " ", initial = "", "No Map Object has been
                    provided. An artificial polygon is used for
                    visualization"))
    map_pseudo(
      object = object, indicator = indicator, panelplot = FALSE,
      MSE = MSE, CV = CV
    )
  } else if (!inherits(x = map_obj, what = "sf")) {
    stop(strwrap(prefix = " ", initial = "",
                 "map_obj is not of class 'sf', 'data.frame'  from the
                 sf package"))
  } else {
    if (length(color) != 2 || !is.vector(color)) {
      stop(strwrap(prefix = " ", initial = "",
                   "col needs to be a vector of length 2 defining the starting,
                   mid and upper color of the map-plot"))
    }

    plot_real(object,
      indicator = indicator,
      MSE = MSE,
      CV = CV,
      map_obj = map_obj,
      map_dom_id = map_dom_id,
      map_tab = map_tab,
      col = color,
      viridis_option = viridis_option, 
      scale_points = scale_points,
      return_data = return_data,
      na.color = na.color, 
      guide = guide,
      save_format = save_format, 
      save_file = save_file
    )
  }
}

map_pseudo <- function(object, indicator, panelplot, MSE, CV) {
  x <- y <- id <- value <- NULL # avoid note due to usage in ggplot
  values <- estimators(
    object = object, indicator = indicator,
    MSE = MSE, CV = CV
  )$ind
  indicator <- colnames(values)[-1]

  tplot <- get_polygone(values = values)

  if (panelplot) {
    ggplot(tplot, aes(x = x, y = y)) +
      geom_polygon(aes(
        group = id,
        fill = value
      )) +
      facet_wrap(~variable,
        ncol = ceiling(sqrt(length(unique(tplot$variable))))
      )
  } else {
    for (ind in indicator) {
      print(ggplot(tplot[tplot$variable == ind, ], aes(x = x, y = y)) +
        ggtitle(paste0(ind)) +
        geom_polygon(aes(
          group = id,
          fill = value
        )))
      cat("Press [enter] to continue")
      line <- readline()
    }
  }
}

plot_real <- function(object,
                      indicator = "all",
                      MSE = FALSE,
                      CV = FALSE,
                      map_obj = NULL,
                      map_dom_id = NULL,
                      map_tab = NULL,
                      col = col,
                      viridis_option = NULL, 
                      scale_points = NULL,
                      na.color = "transparent",
                      return_data = FALSE,
                      guide = NULL,
                      save_file = NULL,
                      save_format = "pdf") {
  if (!is.null(map_obj) && is.null(map_dom_id)) {
    stop("No Domain ID for the map object is given")
  }
  long <- lat <- group <- NULL



  map_data <- estimators(
    object = object, indicator = indicator,
    MSE = MSE, CV = CV
  )$ind

  if (!is.null(map_tab)) {
    names(map_tab)[2] <- map_dom_id
    map_data <- merge(
      x = map_data, y = map_tab,
      by.x = "Domain", by.y = names(map_tab)[1]
    )
    matcher <- match(
      as.vector(unlist(sf::st_drop_geometry(map_obj[map_dom_id]))),
      map_data[, names(map_tab)[2]]
    )

    if (any(is.na(matcher))) {
      if (all(is.na(matcher))) {
        stop(strwrap(prefix = " ", initial = "",
                     "Domains of map_tab and Map object do not match. Check
                     map_tab"))
      } else {
        warning(strwrap(prefix = " ", initial = "",
                         "Not all Domains of map_tab and Map objects could be
                         matched. Check map_tab"))
        #missing_domains <- TRUE
      }
    }
    map_data <- map_data[matcher, ]

    # map_data <- map_data[, !colnames(map_data) %in% c(
    #   map_dom_id,
    #   names(map_tab)
    # ), drop = F]
  } else {
    matcher <- match(as.vector(unlist(sf::st_drop_geometry(map_obj[map_dom_id][, 1]))),
                     map_data[, "Domain"])

    if (any(is.na(matcher))) {
      if (all(is.na(matcher))) {
        stop(strwrap(prefix = " ", initial = "",
                     "Domain of povmap and Map object do not match. Try using
                     map_tab"))
      } else {
        warning(strwrap(prefix = " ", initial = "",
                        "Not all Domains of povmap and Map objects could be matched."))
        #missing_domains <- TRUE 
      }
    }
    map_data <- map_data[matcher, ]
  }

  map_obj <- merge(x = map_obj, y = map_data,by.x=map_dom_id,by.y="Domain",all.x=T)

  indicator <- colnames(map_data)
  indicator <- indicator[!(indicator %in% c("Domain", map_dom_id))]

  for (ind in indicator) {
    map_obj2 <- sf::st_drop_geometry(map_obj)
    map_obj2[ind][,1][!is.finite(map_obj2[ind][,1])] <- NA

    scale_point <- get_scale_points(map_obj2[ind][, 1], ind, scale_points)

    if (is.null(viridis_option)) {
      print(ggplot(data=map_obj,
                   aes(fill = get(ind),color="black")) +
              labs(x = "", y = "", fill=gsub(pattern = "_", replacement = " ", x = ind)) +
              #ggtitle(gsub(pattern = "_", replacement = " ", x = ind)) +
              geom_sf(color = "black") +
              scale_fill_gradient(
                low = col[1], high = col[2],
                limits = scale_point, guide = guide, na.value=na.color
              ) +
              theme(
                axis.ticks = element_blank(), axis.text = element_blank(),
              ) + geom_sf_pattern(data=map_obj[(is.na(map_obj$Head_Count)),], aes(fill=get(ind), colour=""),pattern_color="black") +
              geom_sf(data=map_obj[(is.na(map_obj[,ind])),],color="black",fill="transparent")+
              guides(colour=guide_legend(order=2,"No estimates", override.aes=list(fill=na.color,color="black"))) +
              guides(fill = guide_colourbar(order=1))
      )
    
    }
    else {
      
      print(ggplot(map_obj,
                   aes(fill = get(ind),color="black")) +
              labs(x = "", y = "",fill=gsub(pattern = "_", replacement = " ", x = ind)) +
            
              #ggtitle(gsub(pattern = "_", replacement = " ", x = ind)) +
              geom_sf(color = "black") +
              scale_fill_viridis(option=viridis_option,na.value=na.color) +
              scale_colour_manual(values=NA) +
              theme(
                axis.ticks = element_blank(), axis.text = element_blank()
              ) +
              geom_sf_pattern(data=map_obj[(is.na(map_obj[,ind])),], aes(fill=get(ind), colour=""),pattern_color="black") +
              guides(colour=guide_legend(order=2,"No estimates", override.aes=list(fill=na.color,color="black"))) +
              geom_sf(data=map_obj[(is.na(map_obj[,ind])),],color="black",fill="transparent")+
              guides(fill = guide_colourbar(order=1)) 
      )
      
      
    }
    
    
    if (!is.null(save_file)) {
    ggplot2:::ggsave(file=paste0(save_file,"_",ind,".",save_format))
    }
           
    if (!ind == tail(indicator, 1)) {
      cat("Press [enter] to continue")
      line <- readline()
    }
  }
  if (return_data) {
    return(map_obj)
  }
}

get_polygone <- function(values) {
  if (is.null(dim(values))) {
    values <- as.data.frame(values)
  }
  n <- nrow(values)
  cols <- ceiling(sqrt(n))
  n <- cols^2

  values["id"] <- seq_len(nrow(values))

  poly <- data.frame(
    id = rep(seq_len(n), each = 4),
    ordering = seq_len((n * 4)),
    x = c(0, 1, 1, 0) + rep(0:(cols - 1), each = (cols * 4)),
    y = rep(c(0, 0, 1, 1) + rep(0:(cols - 1), each = 4), cols)
  )

  combo <- merge(poly, values, by = "id", all = TRUE, sort = FALSE)
  melt(combo[order(combo$ordering), ], id.vars = c("id", "x", "y", "ordering"))
}

get_scale_points <- function(y, ind, scale_points) {
  result <- NULL
  if (!is.null(scale_points)) {
    if (is.numeric(scale_points) && length(scale_points) == 2) {
      result <- scale_points
    } else {
      splt <- strsplit(ind, "_\\s*(?=[^_]+$)", perl = TRUE)[[1]]
      indicator_name <- splt[1]
      if (length(splt) == 2) {
        measure <- splt[2]
      } else {
        measure <- "ind"
      }
      if (indicator_name %in% names(scale_points)) {
        pointset <- scale_points[[indicator_name]]
        try(result <- pointset[[measure]])
      }
      if (is.null(result) || length(result) != 2) {
        warning(strwrap(prefix = " ", initial = "",
                        "scale_points is of no apropriate form, default values
                        will be used. See the descriptions and examples for
                        details"))
        result <- NULL
      }
    }
  }
  if (is.null(result)) {
    rg <- range(y, na.rm = TRUE)
    result <- rg
  }
  return(result)
}
