!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE potential_3drism(rismt, vrs, rhog, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... setup potentials for 3D-RISM, which are derived from DFT
  ! ...
  ! ... Variables:
  ! ...   vrs:  DFT's coulomb potential in R-space
  ! ...   rhog: DFT's electronic density in G-space (only for Laue-RISM)
  !
  USE cell_base,      ONLY : alat, tpiba2
  USE constants,      ONLY : pi, e2
  USE control_flags,  ONLY : gamma_only
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE fft_base,       ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect,          ONLY : nl
  USE kinds,          ONLY : DP
  USE lauefft,        ONLY : inv_lauefft_1z, inv_lauefft_2xy
  USE rism,           ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvmol,        ONLY : get_nuniq_in_solVs, solVs, &
                           & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: vrs(1:*)
  COMPLEX(DP),     INTENT(IN)    :: rhog(1:*)
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: nq
  INTEGER                  :: iq
  INTEGER                  :: iiq
  INTEGER                  :: iv
  INTEGER                  :: isolV
  INTEGER                  :: iatom
  REAL(DP)                 :: qv
  REAL(DP),    ALLOCATABLE :: vrss(:)    ! potential(V) in R-Space for Smooth-FFT
  REAL(DP),    ALLOCATABLE :: vrss_s(:)  ! potential(V) in R-Space for Smooth-FFT (Short-range)
  REAL(DP),    ALLOCATABLE :: vrss_l(:)  ! potential(V) in R-Space for Smooth-FFT (Long-range)
  COMPLEX(DP), ALLOCATABLE :: vgss_l(:)  ! potential(V) in G-Space for Smooth-FFT (Long-range)
  COMPLEX(DP), ALLOCATABLE :: vlss_l(:)  ! potential(V) in Laue-rep. for Smooth-FFT (Long-range)
  COMPLEX(DP), ALLOCATABLE :: vright(:)  ! potential coeff. of right-side (for Laue-RISM)
  COMPLEX(DP), ALLOCATABLE :: vleft(:)   ! potential coeff. of left-side  (for Laue-RISM)
  COMPLEX(DP), ALLOCATABLE :: rhogss(:)  ! density(Rho) in G-Space for Smooth-FFT
  COMPLEX(DP), ALLOCATABLE :: aux(:)     ! AUXiliary data for Dense-FFT
  COMPLEX(DP), ALLOCATABLE :: auxs(:)    ! AUXiliary data for Smooth-FFT
#if defined (__DEBUG_RISM)
  !
  CALL start_clock('3DRISM_dft')
#endif
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ng < rismt%cfft%ngmt) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    IF (rismt%ngxy < rismt%lfft%ngxy) THEN
      ierr = IERR_RISM_INCORRECT_DATA_TYPE
      RETURN
    END IF
  END IF
  !
  ! ...
  ! ... allocate working memory
  ! ...
  CALL allocate_works()
  !
  ! ...
  ! ... calculate potential
  ! ...
  ! ... convert coulomb potential: DFT -> 3D-RISM
  CALL interpolate_potential()
  !
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    !
    ! ... convert electronic density: DFT -> 3D-RISM
    CALL interpolate_density()
    !
    ! ... calculate potential of ESM(BC1)
    CALL potential_esm(ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 1
    END IF
  END IF
  !
  ! ...
  ! ... set coulomb potential
  ! ...
  ! ... short-range (R-space)
  IF (rismt%cfft%dfftt%nnr > 0) THEN
    rismt%vsr(:) = 0.0_DP
    rismt%vsr(1:rismt%cfft%dfftt%nnr) = -vrss_s(1:rismt%cfft%dfftt%nnr)
  END IF
  !
  ! ... long-range (R-space)
  IF (rismt%cfft%dfftt%nnr > 0) THEN
    rismt%vlr(:) = 0.0_DP
    rismt%vlr(1:rismt%cfft%dfftt%nnr) = -vrss_l(1:rismt%cfft%dfftt%nnr)
  END IF
  !
  ! ... long-range (G-space)
  IF (rismt%itype == ITYPE_3DRISM) THEN
    IF (rismt%cfft%ngmt > 0) THEN
      rismt%vlgz(:) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      rismt%vlgz(1:rismt%cfft%ngmt) = -vgss_l(1:rismt%cfft%ngmt)
    END IF
  END IF
  !
  ! ... long-range (Laue-rep.)
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    ! ... inside of cell
    IF ((rismt%nrzl * rismt%lfft%ngxy) > 0) THEN
      rismt%vlgz(:) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      rismt%vlgz(1:(rismt%nrzl * rismt%lfft%ngxy)) = -vlss_l(1:(rismt%nrzl * rismt%lfft%ngxy))
    END IF
    !
    ! ... outside of cell
    IF (rismt%lfft%ngxy > 0) THEN
      rismt%vright(:) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      rismt%vleft( :) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      rismt%vright(1:rismt%lfft%ngxy) = -vright(1:rismt%lfft%ngxy)
      rismt%vleft( 1:rismt%lfft%ngxy) = -vleft( 1:rismt%lfft%ngxy)
    END IF
    !
    CALL check_esm_outside(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 1
    END IF
  END IF
  !
  ! ...
  ! ... calculate potential for each solvent's site
  ! ...
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    qv    = solVs(isolV)%charge(iatom)
    !
    ! ... Lennard-Jones potential and short-range coulomb potential (R-space)
    IF (rismt%cfft%dfftt%nnr > 0) THEN
      rismt%usr(:, iiq) = 0.0_DP
      rismt%usr(1:rismt%cfft%dfftt%nnr, iiq) = &
      & rismt%uljr(1:rismt%cfft%dfftt%nnr, iiq) + qv * rismt%vsr(1:rismt%cfft%dfftt%nnr)
    END IF
    !
    ! ... add repulsive-wall potential (R-space)
    IF (rismt%itype == ITYPE_LAUERISM) THEN
      rismt%usr(1:rismt%cfft%dfftt%nnr, iiq) = &
      & rismt%usr(1:rismt%cfft%dfftt%nnr, iiq) + rismt%uwr(1:rismt%cfft%dfftt%nnr, iiq)
    END IF
    !
    IF (rismt%itype == ITYPE_3DRISM) THEN
      !
      ! ... long-range coulomb potential (R-space)
      IF (rismt%cfft%dfftt%nnr > 0) THEN
        rismt%ulr(:, iiq) = 0.0_DP
        rismt%ulr(1:rismt%cfft%dfftt%nnr, iiq) = qv * rismt%vlr(1:rismt%cfft%dfftt%nnr)
      END IF
      !
      ! ... long-range coulomb potential (G-space)
      IF (rismt%cfft%ngmt > 0) THEN
        rismt%ulgz(:, iiq) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
        rismt%ulgz(1:rismt%cfft%ngmt, iiq) = qv * rismt%vlgz(1:rismt%cfft%ngmt)
      END IF
      !
    END IF
    !
  END DO
  !
  ! ... short-range potential (Gxy = 0)
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    IF (rismt%nrzl * rismt%nsite > 0) THEN
      rismt%usg0 = 0.0_DP
    END IF
    !
    CALL corrgxy0_laue(rismt, .TRUE., rismt%usr, rismt%usg0, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 1
    END IF
  END IF
  !
  ! ...
  ! ... normally done
  ! ...
  ierr = IERR_RISM_NULL
  !
1 CONTINUE
  !
  ! ...
  ! ... deallocate working memory
  ! ...
  CALL deallocate_works()
#if defined (__DEBUG_RISM)
  !
  CALL stop_clock('3DRISM_dft')
#endif
  !
CONTAINS
  !
  SUBROUTINE allocate_works()
    !
    IMPLICIT NONE
    !
    ! ... potential
    IF (rismt%cfft%dfftt%nnr > 0) THEN
      ALLOCATE(vrss  (rismt%cfft%dfftt%nnr))
      ALLOCATE(vrss_s(rismt%cfft%dfftt%nnr))
      ALLOCATE(vrss_l(rismt%cfft%dfftt%nnr))
    END IF
    IF (rismt%cfft%ngmt > 0) THEN
      ALLOCATE(vgss_l(rismt%cfft%ngmt))
    END IF
    IF (rismt%itype == ITYPE_LAUERISM) THEN
      IF ((rismt%nrzl * rismt%lfft%ngxy) > 0) THEN
        ALLOCATE(vlss_l(rismt%nrzl * rismt%lfft%ngxy))
      END IF
      IF (rismt%lfft%ngxy > 0) THEN
        ALLOCATE(vright(rismt%lfft%ngxy))
        ALLOCATE(vleft( rismt%lfft%ngxy))
      END IF
    END IF
    !
    ! ... density
    IF (rismt%itype == ITYPE_LAUERISM) THEN
      IF (rismt%cfft%ngmt > 0) THEN
        ALLOCATE(rhogss(rismt%cfft%ngmt))
      END IF
    END IF
    !
    ! ... auxiliary
    IF (dfftp%nnr > 0) THEN
      ALLOCATE(aux   (dfftp%nnr))
    END IF
    IF (rismt%cfft%dfftt%nnr > 0) THEN
      ALLOCATE(auxs  (rismt%cfft%dfftt%nnr))
    END IF
    !
  END SUBROUTINE allocate_works
  !
  SUBROUTINE deallocate_works()
    !
    IMPLICIT NONE
    !
    ! ... potential
    IF (rismt%cfft%dfftt%nnr > 0) THEN
      DEALLOCATE(vrss)
      DEALLOCATE(vrss_s)
      DEALLOCATE(vrss_l)
    END IF
    IF (rismt%cfft%ngmt > 0) THEN
      DEALLOCATE(vgss_l)
    END IF
    IF (rismt%itype == ITYPE_LAUERISM) THEN
      IF ((rismt%nrzl * rismt%lfft%ngxy) > 0) THEN
        DEALLOCATE(vlss_l)
      END IF
      IF (rismt%lfft%ngxy > 0) THEN
        DEALLOCATE(vright)
        DEALLOCATE(vleft)
      END IF
    END IF
    !
    ! ... density
    IF (rismt%itype == ITYPE_LAUERISM) THEN
      IF (rismt%cfft%ngmt > 0) THEN
        DEALLOCATE(rhogss)
      END IF
    END IF
    !
    ! ... auxiliary
    IF (dfftp%nnr > 0) THEN
      DEALLOCATE(aux)
    END IF
    IF (rismt%cfft%dfftt%nnr > 0) THEN
      DEALLOCATE(auxs)
    END IF
    !
  END SUBROUTINE deallocate_works
  !
  SUBROUTINE interpolate_potential()
    !
    ! ... interpolate vrs -> vrss, vrss_s, vrss_l, vgss_l
    !
    IMPLICIT NONE
    !
    INTEGER  :: ir
    INTEGER  :: ig
    REAL(DP) :: gg0
    REAL(DP) :: tt0
    REAL(DP) :: exp0
    REAL(DP) :: v0
    !
    ! ... tau * tau
    tt0 = rismt%tau * rismt%tau
    !
    ! ... vrs -> aux
!$omp parallel do default(shared) private(ir)
    DO ir = 1, dfftp%nnr
      aux(ir) = CMPLX(vrs(ir), 0.0_DP, kind=DP)
    END DO
!$omp end parallel do
    IF (dfftp%nnr > 0) THEN
      CALL fwfft('Dense', aux, dfftp)
    END IF
    !
    ! ... aux -> auxs
    IF (rismt%cfft%dfftt%nnr > 0) THEN
      auxs = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    END IF
!$omp parallel do default(shared) private(ig)
    DO ig = 1, rismt%cfft%ngmt
      auxs(rismt%cfft%nlt(ig)) = aux(nl(ig))
    END DO
!$omp end parallel do
    IF (gamma_only) THEN
!$omp parallel do default(shared) private(ig)
      DO ig = rismt%cfft%gstart_t, rismt%cfft%ngmt
        auxs(rismt%cfft%nltm(ig)) = CONJG(auxs(rismt%cfft%nlt(ig)))
      END DO
!$omp end parallel do
    END IF
    !
    ! ... auxs -> vrss
    IF (rismt%cfft%dfftt%nnr > 0) THEN
      CALL invfft('Custom', auxs, rismt%cfft%dfftt)
    END IF
!$omp parallel do default(shared) private(ir)
    DO ir = 1, rismt%cfft%dfftt%nnr
      vrss(ir) = DBLE(auxs(ir))
    END DO
!$omp end parallel do
    !
    ! ... aux -> auxs, vgss_l
    IF (rismt%cfft%dfftt%nnr > 0) THEN
      auxs = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    END IF
!$omp parallel do default(shared) private(ig, gg0, exp0)
    DO ig = 1, rismt%cfft%ngmt
      gg0  = rismt%cfft%ggt(ig) * tpiba2
      exp0 = EXP(-0.25_DP * gg0 * tt0)
      auxs(rismt%cfft%nlt(ig)) = exp0 * aux(nl(ig))
      vgss_l(ig) = auxs(rismt%cfft%nlt(ig))
    END DO
!$omp end parallel do
    IF (gamma_only) THEN
!$omp parallel do default(shared) private(ig)
      DO ig = rismt%cfft%gstart_t, rismt%cfft%ngmt
        auxs(rismt%cfft%nltm(ig)) = CONJG(auxs(rismt%cfft%nlt(ig)))
      END DO
!$omp end parallel do
    END IF
    !
    ! ... modify auxs, vgss_l
    IF (rismt%itype == ITYPE_LAUERISM) THEN
      IF (rismt%cfft%gstart_t > 1) THEN
        v0 = e2 * pi * tt0 * rhog(1)
        vgss_l(1) = vgss_l(1) - v0
        auxs(rismt%cfft%nlt(1)) = auxs(rismt%cfft%nlt(1)) - v0
      END IF
    END IF
    !
    ! ... auxs -> vrss_s, vrss_l
    IF (rismt%cfft%dfftt%nnr > 0) THEN
      CALL invfft('Custom', auxs, rismt%cfft%dfftt)
    END IF
!$omp parallel do default(shared) private(ir)
    DO ir = 1, rismt%cfft%dfftt%nnr
      vrss_s(ir) = vrss(ir) - DBLE(auxs(ir))
      vrss_l(ir) = DBLE(auxs(ir))
    END DO
!$omp end parallel do
    !
  END SUBROUTINE interpolate_potential
  !
  SUBROUTINE interpolate_density()
    !
    ! ... interpolate rhog -> rhogss
    !
    IMPLICIT NONE
    !
    INTEGER  :: ig
    !
!$omp parallel do default(shared) private(ig)
    DO ig = 1, rismt%cfft%ngmt
      rhogss(ig) = rhog(ig)
    END DO
!$omp end parallel do
    !
  END SUBROUTINE interpolate_density
  !
  SUBROUTINE potential_esm(ierr)
    !
    ! ... calculate coulomb potential of ESM(BC1)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(OUT) :: ierr
    !
    INTEGER  :: ig
    !
    ! ... initialize potentials
    IF ((rismt%nrzl * rismt%lfft%ngxy) > 0) THEN
      vlss_l = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    END IF
    IF (rismt%lfft%ngxy > 0) THEN
      vright = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      vleft  = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    END IF
    !
    ! ... vgss_l -> vlss_l
    IF ((rismt%nrzl * rismt%lfft%ngxy) > 0) THEN
      IF (rismt%cfft%dfftt%nnr > 0) THEN
        auxs = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
!$omp parallel do default(shared) private(ig)
      DO ig = 1, rismt%cfft%ngmt
        auxs(rismt%cfft%nlt(ig)) = vgss_l(ig)
      END DO
!$omp end parallel do
      IF (gamma_only) THEN
!$omp parallel do default(shared) private(ig)
        DO ig = 1, rismt%cfft%ngmt
          auxs(rismt%cfft%nltm(ig)) = CONJG(vgss_l(ig))
        END DO
!$omp end parallel do
      END IF
      !
      IF (rismt%cfft%dfftt%nnr > 0) THEN
        CALL inv_lauefft_1z(rismt%lfft, auxs, vlss_l, rismt%nrzl, rismt%lfft%izcell_start)
      END IF
    END IF
    !
    ! ... add hartree potential to vlss_l
    IF ((rismt%nrzl * rismt%lfft%ngxy) > 0) THEN
      CALL potential_esm_hartree(rismt, rhogss, vlss_l, vright, vleft, ierr)
    ELSE
      ierr = IERR_RISM_NULL
    END IF
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
    !
    ! ... add local potential to vlss_l
    IF ((rismt%nrzl * rismt%lfft%ngxy) > 0) THEN
      CALL potential_esm_local(rismt, 1.0_DP / alat, vlss_l, vright, vleft, ierr)
    ELSE
      ierr = IERR_RISM_NULL
    END IF
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
    !
    ! ... vlss_l -> vrss_l
    IF ((rismt%nrzl * rismt%lfft%ngxy) > 0) THEN
      CALL inv_lauefft_2xy(rismt%lfft, vlss_l, rismt%nrzl, rismt%lfft%izcell_start, vrss_l)
    END IF
    !
    ! ... normally done
    ierr = IERR_RISM_NULL
    !
  END SUBROUTINE potential_esm
  !
END SUBROUTINE potential_3drism
