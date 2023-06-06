! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_vadv
!
!> \brief MPAS ocean vertical advection
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the routine for computing
!>  tendencies for vertical advection.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_vadv

   use mpas_timer
   use ocn_constants
   use ocn_config
   use ocn_mesh

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vel_vadv_tend, &
             ocn_vel_vadv_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      velVadvOff  ! on/off switch for vertical advection


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_vadv_tend
!
!> \brief   Computes tendency term for vertical advection
!> \author  Mark Petersen, Darren Engwirda
!> \date    September 2011
!> \details
!>  This routine computes the vertical advection tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vadv_tend(normalVelocity, layerThickEdgeFlux, &
                                vertAleTransportTop, tend, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,     &!< [in] Horizontal velocity
         layerThickEdgeFlux,     &!< [in] Layer thickness at edge
         vertAleTransportTop  !< [in] Vertical velocity on top layer

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::       &
         iEdge, k,     &! loop indices for edge, vertical loops
         iCell, i,     &! loop indices for dual cells
         kmin, kmax,   &! shallowest and deepest active layers
         vert1, vert2   ! neighbor vertex indices on edge

      real (kind=RKIND) :: &
         r_tmp

      real (kind=RKIND), dimension(:), allocatable :: &
         dudzTop, &! du/dz at layer interface
         wAvgTop   ! w at layer interface

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Set return error code and return early if not turned on
      !*** Start relevant timer

      err = 0
      if (velVadvOff) return
      call mpas_timer_start("vel vadv")


      allocate(dudzTop(nVertLevels+1),wAvgTop(nVertLevels+1))
      !$acc enter data create(dudzTop, wAvgTop)

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(verticesOnEdge, cellsOnVertex, &
      !$acc            minLevelEdgeBot, maxLevelEdgeTop, &
      !$acc            minLevelCell, maxLevelCell, &
      !$acc            kiteAreasOnVertex, areaTriangle, &
      !$acc            vertAleTransportTop, normalVelocity, &
      !$acc            layerThickEdgeFlux, edgeMask) &
      !$acc    private(vert1, vert2, k, i, iCell, kmin, kmax, r_tmp)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(vert1, vert2, k, i, iCell, kmin, kmax, r_tmp)
#endif
      do iEdge = 1, nEdgesOwned
         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)

         ! compute du/dz terms at edge locations
         dudzTop(:) = 0.0_RKIND
         do k = kmin+1, kmax
            dudzTop(k) = (normalVelocity(k-1,iEdge) - &
                          normalVelocity(k-0,iEdge) ) / &
            (0.5_RKIND * (layerThickEdgeFlux(k-1,iEdge) + &
                          layerThickEdgeFlux(k-0,iEdge)))
         end do
         ! not zero on edge bottom!!
         dudzTop(kmax+1) = 1.0_RKIND * dudzTop(kmax)

         ! compute sdot velocities on dual cells
         wAvgTop(:) = 0.0_RKIND
         vert1 = verticesOnEdge(1,iEdge)
         vert2 = verticesOnEdge(2,iEdge)
         do i = 1,vertexDegree
            ! remap to edge from vertex 1
            iCell = cellsOnVertex(i,vert1)
            r_tmp = 0.5_RKIND * &
               kiteAreasOnVertex(i,vert1) / areaTriangle(vert1)
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               wAvgTop(k) = wAvgTop(k) + &
                  r_tmp * vertAleTransportTop(k,iCell)
            end do
            ! remap to edge from vertex 2
            iCell = cellsOnVertex(i,vert2)
            r_tmp = 0.5_RKIND * &
               kiteAreasOnVertex(i,vert2) / areaTriangle(vert2)
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               wAvgTop(k) = wAvgTop(k) + &
                  r_tmp * vertAleTransportTop(k,iCell)
            end do
         end do
         
         ! average w*du/dz to layer midplane
         do k = kmin, kmax
            tend(k,iEdge) = &
            tend(k,iEdge) - edgeMask(k,iEdge) * ( &
               0.50_RKIND * wAvgTop(k+0) * dudzTop(k+0) + &
               0.50_RKIND * wAvgTop(k+1) * dudzTop(k+1) )
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      !$acc exit data delete(dudzTop, wAvgTop)
      deallocate(dudzTop,wAvgTop)

      call mpas_timer_stop("vel vadv")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vadv_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_vadv_init
!
!> \brief   Initializes ocean momentum vertical advection
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical velocity advection in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vadv_init(err)!{{{

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** initialize error code and default values

      err = 0
      velVadvOff = .true.

      !*** Set values based on input configuration

      if (config_vert_coord_movement == 'impermeable_interfaces') then
         velVadvOff = .true.
      else
         velVadvOff = .false.
      end if

      if (config_disable_vel_vadv ) velVadvOff = .true.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vadv_init!}}}

!***********************************************************************

end module ocn_vel_vadv

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
