get_all_factors <- function(n)
{
  prime_factor_tables <- lapply(n, function(i)
  {
    if(i == 1) return(table(1L))
    table(as.numeric(gmp::factorize(i)))
  })
  lapply(prime_factor_tables, function(pft)
  {
    powers <- mapply(
      function(name, value) as.numeric(name) ^ seq.int(0L, value),
      names(pft), 
      pft
    )
    power_grid <- expand.grid(powers)
    sort(apply(power_grid, 1, prod)) 
  })
}

