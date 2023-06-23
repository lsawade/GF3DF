module output

  private
  public :: print_source, print_GF

contains

  subroutine print_GF(GF)

    use ctypes, only: t_GF
    type(t_GF) :: GF

    ! Print format parameters
    CHARACTER(LEN=30) :: integerformat = "(A25 I20)"
    CHARACTER(LEN=30) :: realformat    = "(A25 F26.5)"
    CHARACTER(LEN=30) :: expformat     = "(A25 ES30.5)"

    write (*,*) "**************************************************************"
    write (*,*) "************************* HEADER *****************************"
    write (*,*) "**************************************************************"
    write (*,integerformat)   "do_adjacency_search:", GF%do_adjacency_search
    write (*,integerformat)   "nspec:",               GF%nspec
    write (*,integerformat)   "ngllx:",               GF%ngllx
    write (*,integerformat)   "nglly:",               GF%nglly
    write (*,integerformat)   "ngllz:",               GF%ngllz
    write (*,integerformat)   "midx:",                GF%midx
    write (*,integerformat)   "midy:",                GF%midy
    write (*,integerformat)   "midz:",                GF%midz
    write (*,integerformat)   "nsteps:",              GF%nsteps
    write (*,realformat)      "dt:",                  GF%dt
    write (*,realformat)      "hdur:",                GF%hdur
    write (*,realformat)      "tc:",                  GF%tc
    write (*,expformat)       "factor:",              GF%factor
    write (*,integerformat)   "topography:",          GF%topography

    if (GF%topography == 1) then
      write (*,realformat)      "resolution_topo_file:", GF%resolution_topo_file
      write (*,integerformat)   "nx_bathy:",             GF%nx_bathy
      write (*,integerformat)   "ny_bathy:",             GF%ny_bathy
    endif

    write (*,integerformat)   "ellipticity:",         GF%ellipticity

    write (*,*)
    write (*,*) "----------- Array info ----------------"
    write (*,*)
    write(*, integerformat)  "Number of elements:",   size(GF%ibool,4)
    write(*, integerformat)  "Number of GLL:",        size(GF%xyz,  1)
    write(*, integerformat)  "Ellipticity splines #", size(GF%rspl)
    write(*, integerformat)  "Number of stations",    size(GF%displacement, 1)
    write(*, integerformat)  "Number of components",  size(GF%displacement, 2)
    write(*,*)
    write(*,expformat)      "Displacment Min:", minval(GF%displacement)
    write(*,expformat)      "            Max:", maxval(GF%displacement)
    write(*,expformat)      "           Mean:", sum(GF%displacement)/size(GF%displacement)

    write (*,*) "**************************************************************"


  end subroutine print_GF



  subroutine print_source(source, which)

    use ctypes, only: t_source
    type(t_source) :: source
    integer :: which ! 1 cmtsolution, 2 meshinfo, 3 both

    ! Print format parameters
    CHARACTER(LEN=30) :: charformat    = "(A15 A30)"
    CHARACTER(LEN=30) :: integerformat = "(A25 I20)"
    CHARACTER(LEN=30) :: realformat    = "(A14 F10.5)"
    CHARACTER(LEN=30) :: realformat2   = "(A8  F20.8)"
    CHARACTER(LEN=30) :: expformat     = "(A10 ES14.6)"
    CHARACTER(LEN=30) :: expformat2    = "(A8 ES14.6)"
    CHARACTER(LEN=50) :: PDEformat     = "(xA4I4I3I3I3I3F6.2F9.4F10.4F6.1F4.1F4.1xA25)"

    if (source%force .eqv. .true.) then
      write (*,*) "Printing of Force source not yet implemented. Sorry."
    else
      ! Define Source
      if ((which == 1) .or. (which == 3)) then
        write (*,PDEformat) source%pde_desc, source%year, source%month, source%day, &
                              source%hour, source%minute, source%second, &
                              source%pde_lat, source%pde_lon, source%pde_depth, &
                              source%pde_mb, source%pde_ms, source%pde_region
        write (*,charformat) "event name:    ", source%eventname
        write (*,realformat) "time shift:    ", source%time_shift
        write (*,realformat) "half duration: ", source%hdur
        write (*,realformat) "latitude:      ", source%latitude
        write (*,realformat) "longitude:     ", source%longitude
        write (*,realformat) "depth:         ", source%depth
        write (*,expformat)  "Mrr:      ", source%Mrr
        write (*,expformat)  "Mtt:      ", source%Mtt
        write (*,expformat)  "Mpp:      ", source%Mpp
        write (*,expformat)  "Mrt:      ", source%Mrt
        write (*,expformat)  "Mrp:      ", source%Mrp
        write (*,expformat)  "Mtp:      ", source%Mtp

      endif


      if ((which == 2) .or. (which == 3)) then
        write(*,*) "Source Meshinfo:"
        write(*,*) "------------------------------------------------"
        write(*,expformat2)   "Mxx:    ", source%Mxx
        write(*,expformat2)   "Myy:    ", source%Myy
        write(*,expformat2)   "Mzz:    ", source%Mzz
        write(*,expformat2)   "Mxy:    ", source%Mxy
        write(*,expformat2)   "Mxz:    ", source%Mxz
        write(*,expformat2)   "Myz:    ", source%Myz
        write(*,realformat2) "x:      ", source%x
        write(*,realformat2) "y:      ", source%y
        write(*,realformat2) "z:      ", source%z
        write(*,realformat2) "xix:    ", source%xix
        write(*,realformat2) "xiy:    ", source%xiy
        write(*,realformat2) "xiz:    ", source%xiz
        write(*,realformat2) "etax:   ", source%etax
        write(*,realformat2) "etay:   ", source%etay
        write(*,realformat2) "etaz:   ", source%etaz
        write(*,realformat2) "gammax: ", source%gammax
        write(*,realformat2) "gammay: ", source%gammay
        write(*,realformat2) "gammaz: ", source%gammaz
        write(*,*) "------------------------------------------------"
      endif

    endif

  end subroutine print_source


end module output