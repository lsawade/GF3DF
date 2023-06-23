submodule (sources) io

  implicit none

contains

  ! ==========================================================================
  ! ==========================================================================
  ! ==========================================================================

  module subroutine print_source(source, which)

    ! use types, only: t_source

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

  ! ==========================================================================
  ! ==========================================================================
  ! ==========================================================================

  module function read_cmt(filename) result(sources)

    ! Import statements
    ! use sources_types, only: t_source
    use constants, only: IIN,IMAIN,PI,GRAV,MAX_STRING_LEN,R_PLANET,RHOAV, &
                         NLINES_PER_CMTSOLUTION_SOURCE
    use utils, only: is_digit, is_numeric
    use calendar, only: julian_day

    ! In
    character(len=*) :: filename

    ! Local
    integer :: ios,icounter,nline
    integer :: i, isource, itype, istart, iend, ier
    integer :: NSOURCES
    integer :: yr,mo,da,jda,ho,mi
    double precision :: sec

    character(len=256) :: string, dummystring
    character(len=256) :: eventname

    ! Out
    type(t_source), dimension(:), allocatable :: sources

    ! Number of lines per solution
    nline = NLINES_PER_CMTSOLUTION_SOURCE


    ! Open file to read
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ios)

    ! Check whether the number of lines is a multiple of the CMTSOLUTION format
    icounter = 0
    do while(ios == 0)
      read(IIN,"(a)",iostat=ios) dummystring
      if (ios == 0) icounter = icounter + 1
    enddo

    ! Close source file
    close(IIN)

    if (mod(icounter,nline) /= 0) then
      stop 'total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'
    endif

    ! Get number of sources
    NSOURCES = icounter / nline

    if (NSOURCES < 1) then
      print *,'Error: ',trim(filename),' has ', icounter, 'lines but need ', &
             nline, 'per source... NSOURCES: ', NSOURCES
      stop 'need at least one source in CMTSOLUTION or FORCESOLUTION file'
    endif

    ! Just for now
    allocate(sources(NSOURCES))

    ! Open file again for reading sources
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
    if (ier /= 0) stop 'Error opening CMTSOLUTION file (read_cmt)'

    ! read source number isource
    do isource = 1,NSOURCES

      ! gets header line
      read(IIN,"(a256)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading header line in source ',isource
        stop 'Error reading header line in station in CMTSOLUTION file'
      endif

      ! skips empty lines
      do while( len_trim(string) == 0 )
        read(IIN,"(a256)",iostat=ier) string
        if (ier /= 0) then
          write(IMAIN,*) 'Error reading header line in source ',isource
          stop 'Error reading header blank lines in station in CMTSOLUTION file'
        endif
      enddo

      ! debug
      ! print *,'line ----',string,'----'

      ! Reads the very front
      read(string(1:4),*) sources(isource)%pde_desc
      sources(isource)%pde_desc = trim(sources(isource)%pde_desc)
      ! reads header line with event information (assumes fixed format)
      ! old line: read(string,"(a4,i5,i3,i3,i3,i3,f6.2)") datasource,yr,mo,da,ho,mi,sec

      ! reads header line with event information (free format)
      ! gets rid of the first datasource qualifyer string which can have variable length, like:
      ! "PDE 2014 9 3 .."
      ! " PDEQ2014 9 3 .."
      ! " MLI   1971   1   1 .."
      ! note: globalcmt.org solutions might have missing spaces after datasource qualifier
      !
      ! reads in year,month,day,hour,minutes,seconds
      istart = 1
      do itype = 1,6
        ! determines where first number starts
        do i = istart,len_trim(string)
          if (is_numeric(string(i:i))) then
            istart = i
            exit
          endif
        enddo
        if ( istart >= len_trim(string) ) stop 'Error determining datasource length in header line in CMTSOLUTION file'
        if ( istart <= 1 ) stop 'Error determining datasource length in header line in CMTSOLUTION file'

        ! determines end and length of number
        iend = istart
        do i = istart,len_trim(string)
          if (itype /= 6) then
            ! integer values
            if (.not. is_numeric(string(i:i))) then
              iend = i
              exit
            endif
          else
            ! seconds will have a digit number
            ! digit numbers, e.g. 39.60, can contain '.'
            if (.not. is_digit(string(i:i))) then
              iend = i
              exit
            endif
          endif
        enddo
        iend = iend-1
        if ( iend >= len_trim(string) ) stop 'Error determining number length in header line in CMTSOLUTION file'
        if ( iend < istart ) stop 'Error determining number with negative length in header line in CMTSOLUTION file'

!       ! debug
!       !print *,itype,'line ----',string(istart:iend),'----'

        ! reads in event time information
        ! in case of multiple sources, time refers to the first entry only

        select case (itype)
        case (1)
          ! year (as integer value)
          read(string(istart:iend),*) yr
        case (2)
          ! month (as integer value)
          read(string(istart:iend),*) mo
        case (3)
          ! day (as integer value)
          read(string(istart:iend),*) da
        case (4)
          ! hour (as integer value)
          read(string(istart:iend),*) ho
        case (5)
          ! minutes (as integer value)
          read(string(istart:iend),*) mi
        case (6)
          ! seconds (as float value)
          read(string(istart:iend),*) sec
        end select

        ! advances string
        istart = iend + 1
      enddo

      ! checks time information
      if (yr <= 0 .or. yr > 3000) then
        write(IMAIN,*) 'Error reading year: ',yr,' in source ',isource,'is invalid'
        stop 'Error reading year out of header line in CMTSOLUTION file'
      endif
      if (mo < 1 .or. mo > 12) then
        write(IMAIN,*) 'Error reading month: ',mo,' in source ',isource,'is invalid'
        stop 'Error reading month out of header line in CMTSOLUTION file'
      endif
      if (da < 1 .or. da > 31) then
        write(IMAIN,*) 'Error reading day: ',da,' in source ',isource,'is invalid'
        stop 'Error reading day of header line in CMTSOLUTION file'
      endif
      if (ho < 0 .or. ho > 24) then
        write(IMAIN,*) 'Error reading hour: ',ho,' in source ',isource,'is invalid'
        stop 'Error reading hour of header line in CMTSOLUTION file'
      endif
      if (mi < 0 .or. mi > 59) then
        write(IMAIN,*) 'Error reading minute: ',mi,' in source ',isource,'is invalid'
        stop 'Error reading minute of header line in CMTSOLUTION file'
      endif
      if (sec < 0.0 .or. sec >= 60.0) then
        write(IMAIN,*) 'Error reading second: ',sec,' in source ',isource,'is invalid'
        stop 'Error reading second of header line in CMTSOLUTION file'
      endif

      ! Populate the source
      sources(isource)%year   = yr
      sources(isource)%month  = mo
      sources(isource)%day    = da
      sources(isource)%hour   = ho
      sources(isource)%minute = mi
      sources(isource)%second = sec

      ! gets julian day number
      sources(isource)%jda = julian_day(yr,mo,da)

      ! ignore line with event name
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading event name in source ',isource
        stop 'Error reading event name in station in CMTSOLUTION file'
      endif
      read(string(12:len_trim(string)),*) sources(isource)%eventname

      ! read time shift
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading time shift in source ',isource
        stop 'Error reading time shift in station in CMTSOLUTION file'
      endif
      read(string(12:len_trim(string)),*) sources(isource)%time_shift

      ! read half duration
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading half duration in source ',isource
        stop 'Error reading half duration in station in CMTSOLUTION file'
      endif
      read(string(15:len_trim(string)),*) sources(isource)%hdur

      ! read latitude
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading latitude in source ',isource
        stop 'Error reading latitude in station in CMTSOLUTION file'
      endif
      read(string(10:len_trim(string)),*) sources(isource)%latitude

      ! read longitude
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading longitude in source ',isource
        stop 'Error reading longitude in station in CMTSOLUTION file'
      endif
      read(string(11:len_trim(string)),*) sources(isource)%longitude

      ! read depth
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading depth in source ',isource
        stop 'Error reading depth in station in CMTSOLUTION file'
      endif
      read(string(7:len_trim(string)),*) sources(isource)%depth

      ! seismic moment tensor
      ! CMTSOLUTION: components given in dyne-cm
      ! read Mrr
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading Mrr in source ',isource
        stop 'Error reading Mrr in station in CMTSOLUTION file'
      endif
      read(string(5:len_trim(string)),*) sources(isource)%Mrr

      ! read Mtt
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading Mtt in source ',isource
        stop 'Error reading Mtt in station in CMTSOLUTION file'
      endif
      read(string(5:len_trim(string)),*) sources(isource)%Mtt

      ! read Mpp
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading Mpp in source ',isource
        stop 'Error reading Mpp in station in CMTSOLUTION file'
      endif
      read(string(5:len_trim(string)),*)  sources(isource)%Mpp

      ! read Mrt
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading Mrt in source ',isource
        stop 'Error reading Mrt in station in CMTSOLUTION file'
      endif
      read(string(5:len_trim(string)),*)  sources(isource)%Mrt

      ! read Mrp
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading Mrp in source ',isource
        stop 'Error reading Mrp in station in CMTSOLUTION file'
      endif
      read(string(5:len_trim(string)),*)  sources(isource)%Mrp

      ! read Mtp
      read(IIN,"(a)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading Mtp in source ',isource
        stop 'Error reading Mtp in station in CMTSOLUTION file'
      endif
      read(string(5:len_trim(string)),*)  sources(isource)%Mtp

    enddo

    close(IIN)

    ! !
    ! ! scale and non-dimensionalize the moment tensor
    ! ! CMTSOLUTION file values are in dyne.cm
    ! ! 1 dyne is 1 gram * 1 cm / (1 second)^2
    ! ! 1 Newton is 1 kg * 1 m / (1 second)^2
    ! ! thus 1 Newton = 100,000 dynes
    ! ! therefore 1 dyne.cm = 1e-7 Newton.m
    ! !
    !   scaleM = 1.d7 * RHOAV * (R_PLANET**5) * PI*GRAV*RHOAV
    !   moment_tensor(:,:) = moment_tensor(:,:) / scaleM

  end function read_cmt

end submodule io