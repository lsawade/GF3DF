!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

module get_cmt

  implicit none

  private
  public :: &
    get_cmt_scalar_moment, &
    get_cmt_moment_magnitude, &
    get_cmt_moment_magnitude_from_M0

contains

  double precision function get_cmt_scalar_moment(Mxx,Myy,Mzz,Mxy,Mxz,Myz)

  ! calculates scalar moment (M0)

  use constants, only: PI,GRAV, R_PLANET,RHOAV

  implicit none

  double precision, intent(in) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  ! local parameters
  double precision :: scalar_moment,scaleM

  ! scalar moment:
  ! see equation (1.4) in P.G. Silver and T.H. Jordan, 1982,
  ! "Optimal estimation of scalar seismic moment",
  ! Geophys. J.R. astr. Soc., 70, 755 - 787
  !
  ! or see equation (5.91) in Dahlen & Tromp (1998)
  !
  ! moment tensor M is a symmetric 3x3 tensor, and has six independent components
  !
  ! the euclidean matrix norm is invariant under rotation.
  ! thus, input can be:
  !   Mxx,Myy,Mzz,Mxy,Mxz,Myz
  ! or
  !   Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
  !
  ! euclidean (or Frobenius) norm of a matrix: M0**2 = sum( Mij**2 )
  scalar_moment = Mxx**2 + Myy**2 + Mzz**2 + 2.d0 * ( Mxy**2 + Mxz**2 + Myz**2 )

  ! adds 1/2 to be coherent with double couple or point sources
  scalar_moment = dsqrt(0.5d0*scalar_moment)

  ! note: moment tensor is non-dimensionalized
  !
  ! re-adds scale factor for the moment tensor
  ! CMTSOLUTION file values are in dyne.cm
  ! 1 dyne is 1 gram * 1 cm / (1 second)^2
  ! 1 Newton is 1 kg * 1 m / (1 second)^2
  ! thus 1 Newton = 100,000 dynes
  ! therefore 1 dyne.cm = 1e-7 Newton.m
  scaleM = 1.d7 * RHOAV * (R_PLANET**5) * PI * GRAV * RHOAV

  ! return value (in dyne-cm)
  get_cmt_scalar_moment = scalar_moment * scaleM

  end function

!
!-------------------------------------------------------------------------------------------------
!

  double precision function get_cmt_moment_magnitude(Mxx,Myy,Mzz,Mxy,Mxz,Myz)

  ! calculates moment magnitude (Mw)

  implicit none

  double precision, intent(in) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  ! local parameters
  double precision :: M0,Mw

  ! scalar moment
  M0 = get_cmt_scalar_moment(Mxx,Myy,Mzz,Mxy,Mxz,Myz)

  ! moment magnitude
  Mw = get_cmt_moment_magnitude_from_M0(M0)

  ! return value
  get_cmt_moment_magnitude = Mw

  end function


!
!-------------------------------------------------------------------------------------------------
!

  double precision function get_cmt_moment_magnitude_from_M0(M0)

  ! calculates moment magnitude (Mw) from seismic moment M0

  implicit none

  double precision, intent(in) :: M0

  ! local parameters
  double precision :: Mw

  ! moment magnitude by Hanks & Kanamori, 1979
  ! Mw = 2/3 log( M0 ) - 10.7       (dyne-cm)
  !
  ! alternative forms:
  ! Mw = 2/3 ( log( M0 ) - 16.1 )   (N-m) "moment magnitude" by Hanks & Kanamori(1979) or "energy magnitude" by Kanamori (1977)
  !
  ! Aki & Richards ("Quantitative Seismology",2002):
  ! Mw = 2/3 ( log( M0 ) - 9.1 )    (N-m)
  !
  ! conversion: dyne-cm = 10**-7 N-m
  !
  ! we follow here the USGS magnitude policy:
  ! "Another source of confusion is the form of the formula for converting from scalar moment M0 to moment magnitude, M.
  !  The preferred practice is to use M = (log Mo)/1.5-10.7, where Mo is in dyne-cm (dyne-cm=10-7 N-m),
  !  the definition given by Hanks and Kanamori in 1979. An alternate form in Hanks and Kanamori's paper, M=(log M0-16.1)/1.5,
  !  is sometimes used, with resulting confusion. These formulae look as if they should yield the same result, but the latter
  !  is equivalent to M = (log Mo)/1.5-10.7333. The resulting round-off error occasionally leads to differences of 0.1
  !  in the estimates of moment magnitude released by different groups.
  !  All USGS statements of moment magnitude should use M = (log Mo)/1.5 - 10.7 = 2/3 (log M0) - 10.7 (Hanks & Kanamori, 1979)
  !  for converting from scalar moment M0 to moment magnitude. (..)"
  ! see: http://earthquake.usgs.gov/aboutus/docs/020204mag_policy.php
  !      https://web.archive.org/web/20160428095841/http://earthquake.usgs.gov:80/aboutus/docs/020204mag_policy.php

  if (M0 > 0.0d0) then
    ! this is to ensure M0>0.0 inorder to avoid arithmetic error.
    Mw = 2.d0/3.d0 * log10( max(M0,tiny(M0)) ) - 10.7
  else
    Mw = 0.0d0
  endif

  ! return value
  get_cmt_moment_magnitude_from_M0 = Mw

  end function


end module get_cmt