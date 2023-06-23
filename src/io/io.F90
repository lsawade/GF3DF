
module io

  use sac, only: write_output_SAC
  use GF, only: read_GF, read_GF
  use cmt, only: read_cmt
  use output, only: print_source, print_GF

  private
  public :: read_GF, print_GF, read_GF, write_output_SAC, read_cmt, print_source

end module io