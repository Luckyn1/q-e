!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE do_lauerism(rismt, maxiter, rmsconv, nbox, eta, charge, lboth, iref, title, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... perform Laue-RISM method
  ! ...
  ! ... Variables:
  ! ...   rismt:   data structure
  ! ...   maxiter: maximum number of iterations
  ! ...   rmsconv: RMS of residual vectors to check convergence
  ! ...   nbox:    box size of MDIIS
  ! ...   eta:     step radius of MDIIS
  ! ...   charge:  total charge of solvent system
  ! ...   lboth:   both-hands calculation, or not
  ! ...   iref:    reference type of potential
  ! ...   title:   subtitle of calculation
  ! ...   ierr:    status of calculation
  !
  USE check_stop,     ONLY : check_stop_now, stopped_by_user
  USE control_flags,  ONLY : iverbosity
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE, IERR_RISM_NOT_CONVERGED
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE lauefft,        ONLY : fw_lauefft_2xy, inv_lauefft_2xy
  USE mdiis,          ONLY : mdiis_type, allocate_mdiis, deallocate_mdiis, update_by_mdiis, reset_mdiis
  USE mp,             ONLY : mp_sum, mp_bcast
  USE rism,           ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,        ONLY : get_nuniq_in_solVs, get_nsite_in_solVs, nsolV, solVs, &
                           & iuniq_to_nsite, iuniq_to_isite, isite_to_isolV
  !
  IMPLICIT NONE
  !
  TYPE(rism_type),  INTENT(INOUT) :: rismt
  INTEGER,          INTENT(IN)    :: maxiter
  REAL(DP),         INTENT(IN)    :: rmsconv
  INTEGER,          INTENT(IN)    :: nbox
  REAL(DP),         INTENT(IN)    :: eta
  REAL(DP),         INTENT(IN)    :: charge
  LOGICAL,          INTENT(IN)    :: lboth
  INTEGER,          INTENT(IN)    :: iref
  CHARACTER(LEN=*), INTENT(IN)    :: title
  INTEGER,          INTENT(OUT)   :: ierr
  !
  INTEGER                  :: nq
  INTEGER                  :: iter
  INTEGER                  :: ngrid
  INTEGER                  :: nsite
  INTEGER                  :: nr
  INTEGER                  :: nrr
  INTEGER                  :: nrl
  INTEGER                  :: nrtot
  LOGICAL                  :: lconv
  REAL(DP)                 :: rmscurr
  REAL(DP)                 :: rmssave
  INTEGER                  :: rmswarn
  REAL(DP),    ALLOCATABLE :: hr  (:,:)
  REAL(DP),    ALLOCATABLE :: gr  (:,:)
  REAL(DP),    ALLOCATABLE :: usr (:,:)
  REAL(DP),    ALLOCATABLE :: csr (:,:)
  REAL(DP),    ALLOCATABLE :: dcsr(:,:)
  TYPE(mdiis_type)         :: mdiist
  ! if mdiist is an automatic variable,
  ! pointers in mdiis_type may not work well.
  SAVE                     :: mdiist
#if defined (__DEBUG_RISM)
  CHARACTER(LEN=5)         :: str
#endif
  !
  INTEGER,     PARAMETER   :: NPRINT      = 10
  INTEGER,     PARAMETER   :: MDIIS_EXT   = 3
  INTEGER,     PARAMETER   :: RMSWARN_MAX = 16
  REAL(DP),    PARAMETER   :: RMS_SMALL   = 0.95_DP
  REAL(DP),    PARAMETER   :: RMS_LARGE   = 2.00_DP
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzs < rismt%cfft%dfftt%nr3) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  nq = get_nuniq_in_solVs()
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set dimensions of effective R-space
  nr = rismt%cfft%dfftt%nr1 * rismt%cfft%dfftt%nr2 &
   & * rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1)
  !
  IF (rismt%lfft%gxystart > 1) THEN
    nrr = MAX(0, izright_end0 - izcell_end   )
    nrl = MAX(0, izcell_start - izleft_start0)
  ELSE
    nrr = 0
    nrl = 0
  END IF
  !
  nrtot = nr + nrr + nrl
  !
  ! ... allocate memory
  IF (nrtot * rismt%nsite > 0) THEN
    ALLOCATE(hr(  nrtot, rismt%nsite))
    ALLOCATE(gr(  nrtot, rismt%nsite))
    ALLOCATE(usr( nrtot, rismt%nsite))
    ALLOCATE(csr( nrtot, rismt%nsite))
    ALLOCATE(dcsr(nrtot, rismt%nsite))
  END IF
  !
  CALL allocate_mdiis(mdiist, nbox, nr * rismt%nsite, eta, MDIIS_EXT)
  !
  ! ... reset conditions
  lconv       = .FALSE.
  rismt%avail = .FALSE.
  rmssave     = 1.0E+99_DP
  rmswarn     = 0
  !
  ! ... start Laue-RISM iteration
  WRITE(stdout, '()')
  IF (LEN_TRIM(title) > 0) THEN
    WRITE(stdout, '(5X,"Laue-RISM Calculation (",A,")")') TRIM(title)
  ELSE
    WRITE(stdout, '(5X,"Laue-RISM Calculation")')
  END IF
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"convergence threshold    =",1PE10.3)') rmsconv
#if defined (__DEBUG_RISM)
  !
  IF (iverbosity > 0) THEN
    CALL write_rism_type(rismt)
  END IF
#endif
  !
  ! ... Cs(r) + Cd(z) -> Csd(r)
  CALL corrdipole_laue(rismt, .FALSE., ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    GOTO 100
  END IF
  !
  ! ... Laue-RISM eq. of long-range around the expanded cell
  CALL eqn_lauelong(rismt, lboth, ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    GOTO 100
  END IF
  !
  DO iter = 1, maxiter
    !
    ! ... stop by user
    IF (check_stop_now()) THEN
      EXIT
    END IF
    !
    ! ... FFT: Cs(r) -> Cs(gxy,z)
    CALL fft_csr_to_cslaue()
    !
    ! ... Laue-RISM eq.: Cs(gxy,z), Cd(z) -> H(gxy,z)
    CALL eqn_lauerism(rismt, lboth, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... correct or normalize H(gxy,z)
    ! ... to guarantee total charge and stoichiometry of solvent system
    CALL eqn_laueshort(rismt, lboth, .TRUE., ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    CALL normalize_lauerism(rismt, charge, .FALSE., ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... FFT: H(gxy/=0,z) -> H(r)
    CALL fft_hlaue_to_hr()
    !
    ! ... add H(gxy=0,z) -> H(r)
    CALL corrgxy0_laue(rismt, .FALSE., rismt%hr, rismt%hg0, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... Closure: H(r) -> G(r)
    CALL closure(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... barrier G(r)
    CALL barrier_gr()
    !
    ! ... Residual: G(r) -> dCs(r)
    CALL make_dcsr()
    !
    ! ... clean date out of physical range
    CALL clean_out_of_range(ngrid)
    CALL mp_sum(ngrid, rismt%mp_site%intra_sitg_comm)
    !
    ! ... calculate RMS
    nsite = get_nsite_in_solVs()
    !
    IF (rismt%nr * rismt%nsite > 0) THEN
      CALL rms_residual(ngrid * nsite, rismt%nr * rismt%nsite, &
                      & dcsr, rmscurr, rismt%intra_comm)
    ELSE
      CALL rms_residual(ngrid * nsite, rismt%nr * rismt%nsite, &
                      & dcsr_, rmscurr, rismt%intra_comm)
    END IF
    !
    IF (rmscurr < rmsconv) THEN
      lconv = .TRUE.
    END IF
    !
    ! ... write data
    IF (iverbosity > 0 .OR. MOD(iter - 1, NPRINT) == 0 .OR. lconv .OR. iter == maxiter) THEN
      WRITE(stdout, '(5X,"iter. #",I6,"  RMS(g-h-1)=",1PE10.3,"  nbox=",I3)') &
      & iter, rmscurr, mdiist%nbox
      FLUSH(stdout)
    END IF
#if defined (__DEBUG_RISM)
    !
    IF (iverbosity > 0) THEN
      CALL write_rism_type(rismt)
      !
      WRITE(str, '(I5)') iter
      CALL print_solvavg(rismt, 'rism1.#' // TRIM(ADJUSTL(str)), ierr)
      IF (ierr /= IERR_RISM_NULL) THEN
        GOTO 100
      END IF
    END IF
#endif
    !
    ! ... converged ?
    IF (lconv) THEN
      EXIT
    END IF
    !
    ! ... check RMS to reset MDIIS
    IF (rmscurr > (RMS_SMALL * rmssave)) THEN
      rmswarn = rmswarn + 1
    ELSE
      rmswarn = 0
    END IF
    IF (rmswarn >= RMSWARN_MAX) THEN
      CALL reset_mdiis(mdiist, .TRUE.)
      rmssave = rmscurr
      rmswarn = 0
    END IF
    !
    IF (rmscurr > (RMS_LARGE * rmssave)) THEN
      CALL reset_mdiis(mdiist, .TRUE.)
      rmssave = rmscurr
      rmswarn = 0
    ELSE
      rmssave = MIN(rmscurr, rmssave)
    END IF
    !
    ! ... MDIIS: dCs(r) -> Cs(r)
    IF (rismt%nr * rismt%nsite > 0) THEN
      ! ... optimize with Cs(r) + Cd(z)
      CALL update_by_mdiis(mdiist, rismt%csdr, dcsr, rismt%intra_comm)
      rismt%csr = rismt%csdr
    ELSE
      CALL update_by_mdiis(mdiist, csr_, dcsr_, rismt%intra_comm)
    END IF
    !
    ! ... extract dipole part: Cs(r) -> Cs(r), Cd(z)
    ! ... also perform: Cs(r) + Cd(z) -> Csd(r)
    CALL corrdipole_laue(rismt, .TRUE., ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
  ! ... end Laue-RISM iteration
  END DO
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"End of Laue-RISM calculation")')
  WRITE(stdout, '()')
  !
  ! ... iteration has been converged ?
  IF (lconv) THEN
    ! ... write convergence message
    WRITE(stdout, '()')
    WRITE(stdout, '(5X,"convergence has been achieved in ",I5," iterations")') iter
    WRITE(stdout, '()')
    FLUSH(stdout)
    !
    ! ... Laue-RISM eq. of short-range around the expanded cell
    CALL eqn_laueshort(rismt, lboth, .FALSE., ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... correct or normalize H(gxy,z) to guarantee total charge of solvent system
    CALL normalize_lauerism(rismt, charge, .TRUE., ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... calculate chemical potential
    CALL chempot_lauerism(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... calculate solvation's charge density, potential and energy
    CALL solvation_lauerism(rismt, charge, iref, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... print chemical potential
    CALL print_chempot_lauerism(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... set conditions
    ierr = IERR_RISM_NULL
    rismt%avail = .TRUE.
    !
  ELSE
    ! ... write NOT convergence message
    WRITE(stdout, '()')
    IF (stopped_by_user) THEN
      WRITE(stdout, '(5X,"convergence NOT achieved: stopped by user")')
    ELSE
      WRITE(stdout, '(5X,"convergence NOT achieved")')
    END IF
    WRITE(stdout, '()')
    FLUSH(stdout)
    !
    ! ... set conditions
    ierr = IERR_RISM_NOT_CONVERGED
    rismt%avail = .FALSE.
  END IF
#if defined (__DEBUG_RISM)
  !
  IF (ierr == IERR_RISM_NULL) THEN
    IF (iverbosity > 0) THEN
      CALL write_rism_type(rismt)
    END IF
    !
    CALL print_solvavg(rismt, 'rism1.#end', ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
  END IF
#endif
  !
  ! ... deallocate memory
100 CONTINUE
  !
  IF (nrtot * rismt%nsite > 0) THEN
    DEALLOCATE(hr)
    DEALLOCATE(gr)
    DEALLOCATE(usr)
    DEALLOCATE(csr)
    DEALLOCATE(dcsr)
  END IF
  !
  CALL deallocate_mdiis(mdiist)
  !
CONTAINS
  !
  SUBROUTINE fft_csr_to_cslaue()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: jsite
    !
    DO isite = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      jsite = isite - rismt%mp_site%isite_start + 1
      IF (rismt%nrzs * rismt%ngxy > 0) THEN
        rismt%csgz(:, jsite) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      !
      IF (rismt%cfft%dfftt%nnr > 0 .AND. (rismt%nrzs * rismt%ngxy) > 0) THEN
        CALL fw_lauefft_2xy(rismt%lfft, rismt%csr(:, jsite), rismt%csgz(:, jsite), rismt%nrzs, 1)
      END IF
    END DO
    !
  END SUBROUTINE fft_csr_to_cslaue
  !
  SUBROUTINE fft_hlaue_to_hr()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: jsite
    !
    DO isite = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      jsite = isite - rismt%mp_site%isite_start + 1
      IF (rismt%cfft%dfftt%nnr > 0) THEN
        rismt%hr(:, jsite) = 0.0_DP
      END IF
      !
      IF (rismt%cfft%dfftt%nnr > 0 .AND. (rismt%nrzs * rismt%ngxy) > 0) THEN
        CALL inv_lauefft_2xy(rismt%lfft, rismt%hgz(:, jsite), rismt%nrzs, 1, rismt%hr(:, jsite))
      END IF
    END DO
    !
  END SUBROUTINE fft_hlaue_to_hr
  !
  SUBROUTINE clean_out_of_range(ngrid)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: ngrid
    !
    INTEGER :: ir
    INTEGER :: idx
    INTEGER :: idx0
    INTEGER :: i3min
    INTEGER :: i3max
    INTEGER :: i1, i2, i3
    INTEGER :: iz
    INTEGER :: mgrid
    !
    mgrid = 0
    !
    idx0 = rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x &
       & * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    !
    i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
    !
!$omp parallel do default(shared) private(ir, idx, i1, i2, i3, iz) reduction(+:mgrid)
    DO ir = 1, rismt%cfft%dfftt%nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= rismt%cfft%dfftt%nr3) THEN
        IF (rismt%nsite > 0) THEN
          rismt%csr (ir, :) = 0.0_DP
          rismt%csdr(ir, :) = 0.0_DP
          rismt%hr  (ir, :) = 0.0_DP
          rismt%gr  (ir, :) = 0.0_DP
          dcsr      (ir, :) = 0.0_DP
        END IF
        CYCLE
      END IF
      !
      idx = idx - (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x) * i3
      i2  = idx / rismt%cfft%dfftt%nr1x
      IF (i2 >= rismt%cfft%dfftt%nr2) THEN
        IF (rismt%nsite > 0) THEN
          rismt%csr (ir, :) = 0.0_DP
          rismt%csdr(ir, :) = 0.0_DP
          rismt%hr  (ir, :) = 0.0_DP
          rismt%gr  (ir, :) = 0.0_DP
          dcsr      (ir, :) = 0.0_DP
        END IF
        CYCLE
      END IF
      !
      idx = idx - rismt%cfft%dfftt%nr1x * i2
      i1  = idx
      IF (i1 >= rismt%cfft%dfftt%nr1) THEN
        IF (rismt%nsite > 0) THEN
          rismt%csr (ir, :) = 0.0_DP
          rismt%csdr(ir, :) = 0.0_DP
          rismt%hr  (ir, :) = 0.0_DP
          rismt%gr  (ir, :) = 0.0_DP
          dcsr      (ir, :) = 0.0_DP
        END IF
        CYCLE
      END IF
      !
      IF (i3 < (rismt%cfft%dfftt%nr3 - (rismt%cfft%dfftt%nr3 / 2))) THEN
        iz = i3 + (rismt%cfft%dfftt%nr3 / 2)
      ELSE
        iz = i3 - rismt%cfft%dfftt%nr3 + (rismt%cfft%dfftt%nr3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      !
      IF (iz > rismt%lfft%izright_end .OR. iz < rismt%lfft%izleft_start) THEN
        IF (rismt%nsite > 0) THEN
          rismt%csr (ir, :) = 0.0_DP
          rismt%csdr(ir, :) = 0.0_DP
          rismt%hr  (ir, :) = -1.0_DP
          rismt%gr  (ir, :) = 0.0_DP
          dcsr      (ir, :) = 0.0_DP
        END IF
        CYCLE
      END IF
      IF (iz < rismt%lfft%izright_start .AND. iz > rismt%lfft%izleft_end) THEN
        IF (rismt%nsite > 0) THEN
          rismt%csr (ir, :) = 0.0_DP
          rismt%csdr(ir, :) = 0.0_DP
          rismt%hr  (ir, :) = -1.0_DP
          rismt%gr  (ir, :) = 0.0_DP
          dcsr      (ir, :) = 0.0_DP
        END IF
        CYCLE
      END IF
      !
      mgrid = mgrid + 1
      !
    END DO
!$omp end parallel do
    !
    ngrid = mgrid
    !
  END SUBROUTINE clean_out_of_range
  !
  SUBROUTINE barrier_gr()
    IMPLICIT NONE
    INTEGER  :: ir
    INTEGER  :: idx
    INTEGER  :: idx0
    INTEGER  :: i3min
    INTEGER  :: i3max
    INTEGER  :: i1, i2, i3
    INTEGER  :: iz
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
    idx0 = rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x &
       & * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    !
    i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
    !
!$omp parallel do default(shared) private(ir, idx, i1, i2, i3, iz)
    DO ir = 1, rismt%cfft%dfftt%nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= rismt%cfft%dfftt%nr3) THEN
        CYCLE
      END IF
      !
      idx = idx - (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x) * i3
      i2  = idx / rismt%cfft%dfftt%nr1x
      IF (i2 >= rismt%cfft%dfftt%nr2) THEN
        CYCLE
      END IF
      !
      idx = idx - rismt%cfft%dfftt%nr1x * i2
      i1  = idx
      IF (i1 >= rismt%cfft%dfftt%nr1) THEN
        CYCLE
      END IF
      !
      IF (i3 < (rismt%cfft%dfftt%nr3 - (rismt%cfft%dfftt%nr3 / 2))) THEN
        iz = i3 + (rismt%cfft%dfftt%nr3 / 2)
      ELSE
        iz = i3 - rismt%cfft%dfftt%nr3 + (rismt%cfft%dfftt%nr3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      !
      IF (rismt%lfft%izright_start <= iz .AND. iz < rismt%lfft%izright_gedge) THEN
        rismt%gr(ir, :) = 0.0_DP
      END IF
      !
      IF (rismt%lfft%izleft_gedge < iz .AND. iz <= rismt%lfft%izleft_end) THEN
        rismt%gr(ir, :) = 0.0_DP
      END IF
      !
    END DO
!$omp end parallel do
    !
  END SUBROUTINE barrier_gr
  !
  SUBROUTINE make_dcsr()
    IMPLICIT NONE
    INTEGER  :: iq
    INTEGER  :: iiq
    INTEGER  :: iv
    INTEGER  :: nv
    INTEGER  :: isolV
    INTEGER  :: natom
    REAL(DP) :: rhov
    REAL(DP) :: rhov1
    REAL(DP) :: rhov2
    REAL(DP) :: rhovt
    !
    rhovt = 0.0_DP
    DO isolV = 1, nsolV
      natom = solVs(isolV)%natom
      rhov1 = solVs(isolV)%density
      rhov2 = solVs(isolV)%subdensity
      rhov  = 0.5_DP * (rhov1 + rhov2)
      rhovt = rhovt + rhov * DBLE(natom)
    END DO
    !
    IF (rhovt <= 0.0_DP) THEN ! will not be occurred
      CALL errore('do_lauerism', 'rhovt is not positive', 1)
    END IF
    !
    DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq   = iq - rismt%mp_site%isite_start + 1
      iv    = iuniq_to_isite(1, iq)
      nv    = iuniq_to_nsite(iq)
      isolV = isite_to_isolV(iv)
      rhov1 = solVs(isolV)%density
      rhov2 = solVs(isolV)%subdensity
      rhov  = 0.5_DP * (rhov1 + rhov2) * SQRT(DBLE(nv))
      !
      IF (rismt%nr > 0) THEN
        dcsr(:, iiq) = (rhov / rhovt) * (rismt%gr(:, iiq) - rismt%hr(:, iiq) - 1.0_DP)
      END IF
    END DO
    !
  END SUBROUTINE make_dcsr
  !
END SUBROUTINE do_lauerism
