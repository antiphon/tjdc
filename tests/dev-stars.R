# Misc stars related checks


if(0) {
  # assigning to a matrix a long cell-wise values
  r <- slice(x, "band", 1)
  v <- st_coordinates(r)
  r1 <- r |> mutate(values = v[,1]+v[,2])
  r1 |> plot()
  as_tibble(r1)
}


if(0) {




}
