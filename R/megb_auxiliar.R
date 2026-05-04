# Internal helper: domain summary info used by megb_em
data_info <- function(dom_name, pop, smp) {
  in_dom    <- unique(smp[[dom_name]])
  total_dom <- unique(pop[[dom_name]])
  oosamp    <- !total_dom %in% in_dom

  list(
    n_surv  = length(smp[[dom_name]]),
    n_pop   = length(pop[[dom_name]]),
    d_out   = sum(oosamp),
    d_in    = length(in_dom),
    d_total = length(total_dom),
    ni_smp  = table(smp[[dom_name]]),
    ni_pop  = table(pop[[dom_name]]),
    dom_name = dom_name
  )
}
