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
SUBROUTINE eqn_lauelong(rismt, lboth, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve long-range part of Laue-RISM equation, which is defined as
  ! ...
  ! ...                        /+inf
  ! ...   hl1(gxy,z) = -beta * | dz' ul2(gxy,z') * x21(gxy,z'-z)
  ! ...                        /-inf
  ! ...
  !
  USE constants, ONLY : K_BOLTZMANN_RY, eps8, tpi
  USE cell_base, ONLY : alat
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs, iuniq_to_isite, isite_to_isolV, isite_to_iatom, solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: lboth  ! both-hands calculation, or not
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: nq
  INTEGER                  :: iq1, iq2
  INTEGER                  :: iiq1, iiq2
  INTEGER                  :: iv2
  INTEGER                  :: isolV2
  INTEGER                  :: iatom2
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  INTEGER                  :: iglxy
  INTEGER                  :: jglxy
  INTEGER                  :: irz
  INTEGER                  :: iz1, iz2
  INTEGER                  :: nzint
  INTEGER                  :: izint1
  INTEGER                  :: izint2
  INTEGER                  :: izdelt
  INTEGER                  :: nzright
  INTEGER                  :: izright_sta
  INTEGER                  :: izright_end
  INTEGER                  :: nzleft
  INTEGER                  :: izleft_sta
  INTEGER                  :: izleft_end
  REAL(DP)                 :: beta
  REAL(DP)                 :: qv2
  REAL(DP)                 :: z
  REAL(DP)                 :: gxy
  REAL(DP)                 :: ggxy
  REAL(DP)                 :: ggxy_old
  REAL(DP),    ALLOCATABLE :: xz(:)
  REAL(DP),    ALLOCATABLE :: xgz(:)
  REAL(DP),    ALLOCATABLE :: xgzt(:,:)
  REAL(DP),    ALLOCATABLE :: yz(:)
  REAL(DP),    ALLOCATABLE :: ygz(:)
  REAL(DP),    ALLOCATABLE :: ygzt(:,:)
  REAL(DP)                 :: rstep
  COMPLEX(DP)              :: zstep
  COMPLEX(DP), ALLOCATABLE :: x21(:,:)
  COMPLEX(DP), ALLOCATABLE :: x21out(:,:)
  COMPLEX(DP), ALLOCATABLE :: vl2(:)
  COMPLEX(DP), ALLOCATABLE :: hl1(:,:)
  !
  COMPLEX(DP), PARAMETER   :: C_ZERO = CMPLX(0.0_DP, 0.0_DP, kind=DP)
  COMPLEX(DP), PARAMETER   :: C_ONE  = CMPLX(1.0_DP, 0.0_DP, kind=DP)
  !
  EXTERNAL :: zgemv
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngs < rismt%lfft%nglxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... set dz (in a.u.)
  rstep = alat * rismt%lfft%zstep
  zstep = CMPLX(rstep, 0.0_DP, kind=DP)
  !
  ! ... set integral regions as index of long Z-stick (i.e. expanded cell)
  izright_sta = rismt%lfft%izright_start
  izright_end = rismt%lfft%nrz
  izleft_sta  = 1
  izleft_end  = rismt%lfft%izleft_end
  !
  ! ... count integral points along Z
  nzright = MAX(izright_end - izright_sta + 1, 0)
  nzleft  = MAX(izleft_end  - izleft_sta  + 1, 0)
  nzint   = nzright + nzleft
  !
  ! ... allocate working memory
  IF (rismt%nrzl > 0) THEN
    ALLOCATE(xz(rismt%nrzl))
  END IF
  IF (rismt%nrzl * rismt%ngs > 0) THEN
    ALLOCATE(xgz(rismt%nrzl * rismt%ngs))
  END IF
  IF (rismt%nrzl * rismt%ngs * rismt%nsite > 0) THEN
    ALLOCATE(xgzt(rismt%nrzl * rismt%ngs, rismt%nsite))
  END IF
  !
  IF (rismt%nrzl > 0) THEN
    ALLOCATE(yz(rismt%nrzl))
  END IF
  IF (lboth) THEN
    IF (rismt%nrzl * rismt%ngs > 0) THEN
      ALLOCATE(ygz(rismt%nrzl * rismt%ngs))
    END IF
    IF (rismt%nrzl * rismt%ngs * rismt%nsite > 0) THEN
      ALLOCATE(ygzt(rismt%nrzl * rismt%ngs, rismt%nsite))
    END IF
  END IF
  !
  IF (nzint > 0) THEN
    ALLOCATE(x21(nzint, nzint))
  END IF
  IF (rismt%lfft%nrz * nzint > 0) THEN
    ALLOCATE(x21out(rismt%lfft%nrz, nzint))
  END IF
  IF (rismt%lfft%nrz > 0) THEN
    ALLOCATE(vl2(rismt%lfft%nrz))
  END IF
  IF (nzint * rismt%lfft%ngxy > 0) THEN
    ALLOCATE(hl1(nzint, rismt%lfft%ngxy))
  END IF
  !
  !                         ----
  ! ... reduce x21 :   x1 = >    q2 * x21
  !                         ----
  !                           2
  DO iq1 = 1, nq
    ! ... properties of site1
    IF (rismt%mp_site%isite_start <= iq1 .AND. iq1 <= rismt%mp_site%isite_end) THEN
      iiq1 = iq1 - rismt%mp_site%isite_start + 1
    ELSE
      iiq1 = 0
    END IF
    !
    ! x21 -> x1 (the right-hand side)
    CALL reduce_xgs(rismt%xgs(:, :, iq1), rismt%nrzl * rismt%ngs)
    IF (iiq1 > 0) THEN
      IF (rismt%nrzl * rismt%ngs > 0) THEN
        xgzt(:, iiq1) = xgz(:)
      END IF
    END IF
    !
    IF (.NOT. lboth) THEN
      CYCLE
    END IF
    !
    ! x21 -> x1 (the left-hand side)
    CALL reduce_xgs(rismt%ygs(:, :, iq1), rismt%nrzl * rismt%ngs)
    IF (iiq1 > 0) THEN
      IF (rismt%nrzl * rismt%ngs > 0) THEN
        ygzt(:, iiq1) = xgz(:)
      END IF
    END IF
    !
  END DO
  !
  ! ... calculate integral equation:
  !
  !                             /+inf
  !       hl1(gxy,z1) = -beta * | dz2 vl(gxy,z2) * x1(gxy,z2-z1)
  !                             /-inf
  !
  DO iq1 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq1 = iq1 - rismt%mp_site%isite_start + 1
    xgz(:) = xgzt(:, iiq1)
    IF (lboth) THEN
      ygz(:) = ygzt(:, iiq1)
    END IF
    !
    ! ... in case of zleft <= z <= zright
    !
    !                             /zright
    !       hl1(gxy,z1) = -beta * | dz2 vl(gxy,z2) * x1(gxy,z2-z1)
    !                             /zleft
    !
    ! ... solve Laue-RISM equation for each igxy
    ggxy_old = -1.0_DP
    !
    DO igxy = 1, rismt%lfft%ngxy
      jgxy  = rismt%nrzl * (igxy - 1)
      iglxy = rismt%lfft%igtonglxy(igxy)
      jglxy = rismt%nrzl * (iglxy - 1)
      ggxy  = rismt%lfft%glxy(iglxy)
      gxy   = rismt%lfft%gnxy(igxy)
      !
      ! ... x(z2-z1)
      IF (ggxy > (ggxy_old + eps8)) THEN
        ggxy_old = ggxy
        !
        IF (.NOT. lboth) THEN
          ! ... susceptibility matrix of one-hand calculation, which is symmetric
          xz(1:rismt%nrzl) = xgz((1 + jglxy):(rismt%nrzl + jglxy))
          !
          DO iz1 = izleft_sta, izleft_end
            izint1 = iz1 - izleft_sta + 1
            DO iz2 = izleft_sta, iz1
              izint2 = iz2 - izleft_sta + 1
              izdelt = iz1 - iz2 + 1
              x21(izint2, izint1) = CMPLX(xz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
          DO iz1 = izright_sta, izright_end
            izint1 = nzleft + iz1 - izright_sta + 1
            DO iz2 = izright_sta, iz1
              izint2 = nzleft + iz2 - izright_sta + 1
              izdelt = iz1 - iz2 + 1
              x21(izint2, izint1) = CMPLX(xz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
          DO iz1 = izright_sta, izright_end
            izint1 = nzleft + iz1 - izright_sta + 1
            DO iz2 = izleft_sta, izleft_end
              izint2 = iz2 - izleft_sta + 1
              izdelt = iz1 - iz2 + 1
              x21(izint2, izint1) = CMPLX(xz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
          DO izint1 = 1, nzint
            DO izint2 = 1, (izint1 - 1)
              x21(izint1, izint2) = x21(izint2, izint1)
            END DO
          END DO
          !
        ELSE !IF (lboth) THEN
          ! ... susceptibility matrix of both-hands calculation
          xz(1:rismt%nrzl) = xgz((1 + jglxy):(rismt%nrzl + jglxy))
          yz(1:rismt%nrzl) = ygz((1 + jglxy):(rismt%nrzl + jglxy))
          !
          DO iz1 = izleft_sta, izleft_end
            izint1 = iz1 - izleft_sta + 1
            DO iz2 = izleft_sta, izleft_end
              izint2 = iz2 - izleft_sta + 1
              izdelt = ABS(iz1 - iz2) + 1
              x21(izint2, izint1) = CMPLX(yz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
          DO iz1 = izright_sta, izright_end
            izint1 = nzleft + iz1 - izright_sta + 1
            DO iz2 = izright_sta, izright_end
              izint2 = nzleft + iz2 - izright_sta + 1
              izdelt = ABS(iz1 - iz2) + 1
              x21(izint2, izint1) = CMPLX(xz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
          DO iz1 = izleft_sta, izleft_end
            izint1 = iz1 - izleft_sta + 1
            DO iz2 = izright_sta, izright_end
              izint2 = nzleft + iz2 - izright_sta + 1
              izdelt = iz2 - iz1 + 1
              x21(izint2, izint1) = CMPLX(yz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
          DO iz1 = izright_sta, izright_end
            izint1 = nzleft + iz1 - izright_sta + 1
            DO iz2 = izleft_sta, izleft_end
              izint2 = iz2 - izleft_sta + 1
              izdelt = iz1 - iz2 + 1
              x21(izint2, izint1) = CMPLX(xz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
        END IF
        !
      END IF
      !
      ! ... vl(z2)
      DO iz2 = izleft_sta, izleft_end
        izint2 = iz2 - izleft_sta + 1
        vl2(izint2) = rismt%vlgz(iz2 + jgxy)
      END DO
      DO iz2 = izright_sta, izright_end
        izint2 = nzleft + iz2 - izright_sta + 1
        vl2(izint2) = rismt%vlgz(iz2 + jgxy)
      END DO
      !
      ! ... hl(z1)
      IF (nzint > 0) THEN
        CALL zgemv('T', nzint, nzint, zstep, x21, nzint, vl2, 1, C_ZERO, hl1(1, igxy), 1)
      END IF
      !
    END DO
    !
    ! ... in case of z < zleft
    !
    !                             /zleft
    !       hl1(gxy,z1) = -beta * | dz2 vl(gxy,z2) * x1(gxy,z2-z1)
    !                             /-inf
    !
    ! ... solve Laue-RISM equation for each igxy
    IF (rismt%lfft%xleft) THEN
      ggxy_old = -1.0_DP
      !
      DO igxy = 1, rismt%lfft%ngxy
        IF (.NOT. rismt%do_vleft(igxy)) THEN
          CYCLE
        END IF
        !
        jgxy  = rismt%nrzl * (igxy - 1)
        iglxy = rismt%lfft%igtonglxy(igxy)
        jglxy = rismt%nrzl * (iglxy - 1)
        ggxy  = rismt%lfft%glxy(iglxy)
        gxy   = rismt%lfft%gnxy(igxy)
        !
        ! ... x(z2-z1)
        IF (ggxy > (ggxy_old + eps8)) THEN
          ggxy_old = ggxy
          x21out = C_ZERO
          !
          xz(1:rismt%nrzl) = xgz((1 + jglxy):(rismt%nrzl + jglxy))
          IF (.NOT. lboth) THEN
            yz(1:rismt%nrzl) = xgz((1 + jglxy):(rismt%nrzl + jglxy))
          ELSE
            yz(1:rismt%nrzl) = ygz((1 + jglxy):(rismt%nrzl + jglxy))
          END IF
          !
          DO iz1 = izleft_sta, izleft_end
            izint1 = iz1 - izleft_sta + 1
            DO iz2 = (iz1 - rismt%lfft%nrz + 1), 0
              izint2 = iz2 + rismt%lfft%nrz
              izdelt = iz1 - iz2 + 1
              x21out(izint2, izint1) = CMPLX(yz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
          DO iz1 = izright_sta, izright_end
            izint1 = nzleft + iz1 - izright_sta + 1
            DO iz2 = (iz1 - rismt%lfft%nrz + 1), 0
              izint2 = iz2 + rismt%lfft%nrz
              izdelt = iz1 - iz2 + 1
              x21out(izint2, izint1) = CMPLX(xz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
        END IF
        !
        ! ... vl(z2)
        IF (rismt%lfft%gxystart > 1 .AND. igxy == 1) THEN
          ! ... Gxy = 0
          DO iz2 = (-rismt%lfft%nrz + 1), 0
            izint2 = iz2 + rismt%lfft%nrz
            z = rismt%lfft%zleft + rismt%lfft%zoffset + rismt%lfft%zstep * DBLE(iz2 - 1)
            vl2(izint2) = CMPLX(DBLE(rismt%vleft(igxy)) * z &
                           & + AIMAG(rismt%vleft(igxy)), 0.0_DP, kind=DP)
          END DO
          !
        ELSE
          ! ... Gxy /= 0
          DO iz2 = (-rismt%lfft%nrz + 1), 0
            izint2 = iz2 + rismt%lfft%nrz
            z = rismt%lfft%zleft + rismt%lfft%zoffset + rismt%lfft%zstep * DBLE(iz2 - 1)
            vl2(izint2) = EXP( tpi * gxy * (z - rismt%lfft%zleft)) * rismt%vleft(igxy)
          END DO
          !
        END IF
        !
        ! ... hl(z1)
        IF (nzint > 0) THEN
          CALL zgemv('T', rismt%lfft%nrz, nzint, zstep, &
                   & x21out, rismt%lfft%nrz, vl2, 1, C_ONE, hl1(1, igxy), 1)
        END IF
        !
      END DO
    END IF
    !
    ! ... in case of zright < z
    !
    !                             /+inf
    !       hl1(gxy,z1) = -beta * | dz2 vl(gxy,z2) * x1(gxy,z2-z1)
    !                             /zright
    !
    ! ... solve Laue-RISM equation for each igxy
    IF (rismt%lfft%xright) THEN
      ggxy_old = -1.0_DP
      !
      DO igxy = 1, rismt%lfft%ngxy
        IF (.NOT. rismt%do_vright(igxy)) THEN
          CYCLE
        END IF
        !
        jgxy  = rismt%nrzl * (igxy - 1)
        iglxy = rismt%lfft%igtonglxy(igxy)
        jglxy = rismt%nrzl * (iglxy - 1)
        ggxy  = rismt%lfft%glxy(iglxy)
        gxy   = rismt%lfft%gnxy(igxy)
        !
        ! ... x(z2-z1)
        IF (ggxy > (ggxy_old + eps8)) THEN
          ggxy_old = ggxy
          x21out = C_ZERO
          !
          xz(1:rismt%nrzl) = xgz((1 + jglxy):(rismt%nrzl + jglxy))
          IF (.NOT. lboth) THEN
            yz(1:rismt%nrzl) = xgz((1 + jglxy):(rismt%nrzl + jglxy))
          ELSE
            yz(1:rismt%nrzl) = ygz((1 + jglxy):(rismt%nrzl + jglxy))
          END IF
          !
          DO iz1 = izleft_sta, izleft_end
            izint1 = iz1 - izleft_sta + 1
            DO iz2 = (rismt%lfft%nrz + 1), (iz1 + rismt%lfft%nrz - 1)
              izint2 = iz2 - rismt%lfft%nrz
              izdelt = iz2 - iz1 + 1
              x21out(izint2, izint1) = CMPLX(yz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
          DO iz1 = izright_sta, izright_end
            izint1 = nzleft + iz1 - izright_sta + 1
            DO iz2 = (rismt%lfft%nrz + 1), (iz1 + rismt%lfft%nrz - 1)
              izint2 = iz2 - rismt%lfft%nrz
              izdelt = iz2 - iz1 + 1
              x21out(izint2, izint1) = CMPLX(xz(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
          !
        END IF
        !
        ! ... vl(z2)
        IF (rismt%lfft%gxystart > 1 .AND. igxy == 1) THEN
          ! ... Gxy = 0
          DO iz2 = (rismt%lfft%nrz + 1), (2 * rismt%lfft%nrz)
            izint2 = iz2 - rismt%lfft%nrz
            z = rismt%lfft%zleft + rismt%lfft%zoffset + rismt%lfft%zstep * DBLE(iz2 - 1)
            vl2(izint2) = CMPLX(DBLE(rismt%vright(igxy)) * z &
                           & + AIMAG(rismt%vright(igxy)), 0.0_DP, kind=DP)
          END DO
          !
        ELSE
          ! ... Gxy /= 0
          DO iz2 = (rismt%lfft%nrz + 1), (2 * rismt%lfft%nrz)
            izint2 = iz2 - rismt%lfft%nrz
            z = rismt%lfft%zleft + rismt%lfft%zoffset + rismt%lfft%zstep * DBLE(iz2 - 1)
            vl2(izint2) = EXP(-tpi * gxy * (z - rismt%lfft%zright)) * rismt%vright(igxy)
          END DO
          !
        END IF
        !
        ! ... hl(z1)
        IF (nzint > 0) THEN
          CALL zgemv('T', rismt%lfft%nrz, nzint, zstep, &
                   & x21out, rismt%lfft%nrz, vl2, 1, C_ONE, hl1(1, igxy), 1)
        END IF
        !
      END DO
    END IF
    !
    ! ... copy hl1 -> hlgz
    IF ((rismt%nrzl * rismt%ngxy) > 0) THEN
      rismt%hlgz(:, iiq1) = C_ZERO
    END IF
    !
    DO igxy = 1, rismt%lfft%ngxy
      jgxy = rismt%nrzl * (igxy - 1)
      DO iz1 = izleft_sta, izleft_end
        izint1 = iz1 - izleft_sta + 1
        rismt%hlgz(iz1 + jgxy, iiq1) = -beta * hl1(izint1, igxy)
      END DO
      DO iz1 = izright_sta, izright_end
        izint1 = nzleft + iz1 - izright_sta + 1
        rismt%hlgz(iz1 + jgxy, iiq1) = -beta * hl1(izint1, igxy)
      END DO
    END DO
    !
  END DO
  !
  ! ... deallocate working memory
  IF (rismt%nrzl > 0) THEN
    DEALLOCATE(xz)
  END IF
  IF (rismt%nrzl * rismt%ngs > 0) THEN
    DEALLOCATE(xgz)
  END IF
  IF (rismt%nrzl * rismt%ngs * rismt%nsite > 0) THEN
    DEALLOCATE(xgzt)
  END IF
  !
  IF (rismt%nrzl > 0) THEN
    DEALLOCATE(yz)
  END IF
  IF (lboth) THEN
    IF (rismt%nrzl * rismt%ngs > 0) THEN
      DEALLOCATE(ygz)
    END IF
    IF (rismt%nrzl * rismt%ngs * rismt%nsite > 0) THEN
      DEALLOCATE(ygzt)
    END IF
  END IF
  !
  IF (nzint > 0) THEN
    DEALLOCATE(x21)
  END IF
  IF (rismt%lfft%nrz * nzint > 0) THEN
    DEALLOCATE(x21out)
  END IF
  IF (rismt%lfft%nrz > 0) THEN
    DEALLOCATE(vl2)
  END IF
  IF (nzint * rismt%lfft%ngxy > 0) THEN
    DEALLOCATE(hl1)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE reduce_xgs(xgs, ndim)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: xgs(1:ndim, 1:*)
    INTEGER,  INTENT(IN) :: ndim
    !
    REAL(DP) :: xg_int
    !
    IF (rismt%nrzl * rismt%ngs > 0) THEN
      xgz = 0.0_DP
    END IF
    !
    DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      ! ... properties of site2
      iiq2   = iq2 - rismt%mp_site%isite_start + 1
      iv2    = iuniq_to_isite(1, iq2)
      isolV2 = isite_to_isolV(iv2)
      iatom2 = isite_to_iatom(iv2)
      qv2    = solVs(isolV2)%charge(iatom2)
      !
      ! ... calculate x21 -> x1
      IF (rismt%nrzl * rismt%ngs > 0) THEN
        xgz(:) = xgz(:) + qv2 * xgs(1:(rismt%nrzl * rismt%ngs), iiq2)
      END IF
      !
    END DO
    !
    IF (rismt%nrzl * rismt%ngs > 0) THEN
      CALL mp_sum(xgz, rismt%mp_site%inter_sitg_comm)
    END IF
    !
    ! ... correct at G = 0
    IF (rismt%lfft%gxystart > 1) THEN
      xg_int = rstep * xgz(1)
      DO irz = 2, rismt%lfft%nrz
        xg_int = xg_int + 2.0_DP * rstep * xgz(irz)
      END DO
      DO irz = 1, rismt%lfft%nrz
        xgz(irz) = xgz(irz) - xg_int / rstep / DBLE(2 * (rismt%lfft%nrz - 1) + 1)
      END DO
    END IF
    !
  END SUBROUTINE reduce_xgs
  !
END SUBROUTINE eqn_lauelong