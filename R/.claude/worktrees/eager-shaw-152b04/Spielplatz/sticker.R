pacman::p_load(hexSticker)

sticker(

  # image
  "man/figures/poverty_mapping_image.jpg", # https://unsplash.com/photos/7wBFsHWQDlk
  s_x=1, # slightly to right to appear centered
  s_y=1,
  s_width=1.5,
  s_height=1.5,
  white_around_sticker = T,

  # package name
  package="povmap",
  p_size=50,
  p_color = "ivory", # 00030A 010101 #383838
  p_y = 1.6,

  # Output file
  filename="man/figures/logo.png",

  # Background colour
  h_fill = "#F0F0F0", # #F0F0F0


  # Border
  # Grey colours: https://www.w3schools.com/colors/colors_shades.asp
  h_color = "ivory",   # 3F4243 7F2B94 3B2691 4238AF
  # h_size = 1.5,

  dpi = 1000 # otherwise the final fantasy image quality is not good
)
