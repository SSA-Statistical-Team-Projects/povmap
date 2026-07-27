
# Generate emdi object with additional indicators; here via function ebp()
emdi_model <- ebp(
  fixed = eqIncome ~ gender + eqsize + cash +
    self_empl + unempl_ben + age_ben + surv_ben + sick_ben +
    dis_ben + rent + fam_allow + house_allow + cap_inv +
    tax_adj, pop_data = eusilcA_pop,
  pop_domains = "district", smp_data = eusilcA_smp,
  smp_domains = "district", threshold = 11064.82,
  transformation = "box.cox", L = 2, MSE = TRUE, B = 2
)


### change the shapefile to an sf object
# shape_austria_dis <- sf::st_as_sf(shape_austria_dis, crs = 4326)
#
# save(shape_austria_dis, file = "inst/shapes/shape_austria_dis.rda")

#### test that the new shapefile works
load_shapeaustria()


# First find the right order
dom_ord <- match(shape_austria_dis$PB, emdi_model$ind$Domain)

# Create the mapping table based on the order obtained above
map_tab <- data.frame(
  pop_data_id = emdi_model$ind$Domain[dom_ord],
  shape_id = shape_austria_dis$BKZ
)

# Create map plot for mean indicator - point and CV estimates but no MSE
# using the numerical domain identifiers of the shape file

plot_real(
  object = emdi_model, MSE = FALSE, CV = TRUE,
  map_obj = shape_austria_dis, indicator = c("Mean"),
  map_dom_id = "BKZ", map_tab = map_tab,
  col = c("lightblue", "darkblue")
)

map_plot(
  object = emdi_model, MSE = FALSE, CV = TRUE,
  map_obj = shape_austria_dis, indicator = c("Mean"),
  map_dom_id = "BKZ", map_tab = map_tab
)

