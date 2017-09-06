!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE xml_io_rism
  !----------------------------------------------------------------------------
  !
  ! ... this module contains subroutines used to read and write
  ! ... 1D- and 3D-RISM data in XML format
  !
  USE iotk_module
  USE constants,   ONLY : eps8
  USE kinds,       ONLY : DP
  USE lauefft,     ONLY : lauefft_type
  USE mp,          ONLY : mp_rank, mp_size, mp_sum, mp_get, mp_bcast, mp_barrier
  USE xml_io_base, ONLY : check_file_exst
  USE parallel_include
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define constants
#if defined (__1DRISM_BINARY)
  LOGICAL :: l1drism_binary = .TRUE.
#else
  LOGICAL :: l1drism_binary = .FALSE.
#endif
  LOGICAL :: l3drism_binary = .TRUE.
  !
  ! ... public components
  PUBLIC :: l1drism_binary
  PUBLIC :: l3drism_binary
  PUBLIC :: read_1drism_xml
  PUBLIC :: write_1drism_xml
  PUBLIC :: read_3drism_xml
  PUBLIC :: write_3drism_xml
  PUBLIC :: read_lauerism_xml
  PUBLIC :: write_lauerism_xml
  PUBLIC :: read_lauedipole_xml
  PUBLIC :: write_lauedipole_xml
  !
CONTAINS
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_1drism_xml(rism1d_file_base, zvv, name, ngrid, &
                            & nsite, ipp, npp, ionode, intra_group_comm)
    !------------------------------------------------------------------------
    !
    ! ... Writes 1D-RISM's correlation function, one site at a time.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: rism1d_file_base
    REAL(DP),         INTENT(IN) :: zvv(:,:)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER,          INTENT(IN) :: ngrid
    INTEGER,          INTENT(IN) :: nsite
    INTEGER,          INTENT(IN) :: ipp(:)
    INTEGER,          INTENT(IN) :: npp(:)
    LOGICAL,          INTENT(IN) :: ionode
    INTEGER,          INTENT(IN) :: intra_group_comm
    !
    INTEGER                 :: ierr
    INTEGER                 :: isite
    INTEGER                 :: rism1d_unit
    CHARACTER(LEN=256)      :: rism1d_file
    CHARACTER(LEN=10)       :: rism1d_extension
    CHARACTER(iotk_attlenx) :: attr
    INTEGER                 :: io_group
    INTEGER                 :: me_group
    REAL(DP), ALLOCATABLE   :: zvv1(:)
    !
    INTEGER, EXTERNAL       :: find_free_unit
    !
    ! ... get process info.
    me_group = mp_rank(intra_group_comm)
    !
    ! ... decide file name and unit
    rism1d_extension = '.dat'
    IF (.NOT. l1drism_binary) THEN
      rism1d_extension = '.xml'
    END IF
    !
    rism1d_file = TRIM(rism1d_file_base) // TRIM(rism1d_extension)
    rism1d_unit = find_free_unit()
    !
    ! ... open file
    IF (ionode) THEN
      CALL iotk_open_write(rism1d_unit, FILE=rism1d_file,  BINARY=l1drism_binary, IERR=ierr)
      CALL errore('write_rism1d_xml', &
      & 'cannot open ' // TRIM(rism1d_file) // ' file for writing', ierr)
    END IF
    !
    ! ... write header
    IF (ionode) THEN
      CALL iotk_write_begin(rism1d_unit, "_1D-RISM")
      !
      CALL iotk_write_attr(attr, "name", TRIM(name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "ngrid", ngrid)
      CALL iotk_write_attr(attr, "nsite", nsite)
      CALL iotk_write_empty(rism1d_unit, "INFO", attr)
    END IF
    !
    ! ... find the index of the ionode
    io_group = 0
    IF (ionode) THEN
      io_group = me_group
    END IF
    CALL mp_sum(io_group, intra_group_comm)
    !
    ! ... write zvv for each site
    ALLOCATE(zvv1(ngrid))
    !
    DO isite = 1, nsite
#if defined (__MPI)
      CALL MPI_GATHERV(zvv(1, isite), npp(me_group + 1), &
         & MPI_DOUBLE_PRECISION, zvv1, npp, ipp, &
         & MPI_DOUBLE_PRECISION, io_group, intra_group_comm, ierr)
      !
      IF (ierr /= MPI_SUCCESS) THEN
        CALL errore('write_1drism_xml', 'error at MPI_GATHERV', 1)
      END IF
#else
      zvv1(1:ngrid) = zvv(1:ngrid, isite)
#endif
      !
      IF (ionode) THEN
        CALL iotk_write_dat(rism1d_unit, "site" // iotk_index(isite), zvv1)
      END IF
    END DO
    !
    DEALLOCATE(zvv1)
    !
    ! ... close file
    IF (ionode) THEN
      CALL iotk_write_end(rism1d_unit, "_1D-RISM")
      !
      CALL iotk_close_write(rism1d_unit)
    END IF
    !
  END SUBROUTINE write_1drism_xml
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_1drism_xml(rism1d_file_base, zvv, ngrid, &
                           & nsite, ipp, npp, ionode, ionode_id, intra_group_comm)
    !------------------------------------------------------------------------
    !
    ! ... Reads 1D-RISM's correlation function, one site at a time.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN)  :: rism1d_file_base
    REAL(DP),         INTENT(OUT) :: zvv(:,:)
    INTEGER,          INTENT(IN)  :: ngrid
    INTEGER,          INTENT(IN)  :: nsite
    INTEGER,          INTENT(IN)  :: ipp(:)
    INTEGER,          INTENT(IN)  :: npp(:)
    LOGICAL,          INTENT(IN)  :: ionode
    INTEGER,          INTENT(IN)  :: ionode_id
    INTEGER,          INTENT(IN)  :: intra_group_comm
    !
    INTEGER                 :: ierr
    INTEGER                 :: isite
    INTEGER                 :: rism1d_unit
    CHARACTER(LEN=256)      :: rism1d_file
    CHARACTER(iotk_attlenx) :: attr
    INTEGER                 :: ngrid_
    INTEGER                 :: nsite_
    INTEGER                 :: io_group
    INTEGER                 :: me_group
    LOGICAL                 :: exst
    REAL(DP), ALLOCATABLE   :: zvv1(:)
    !
    INTEGER, EXTERNAL       :: find_free_unit
    !
    ! ... get process info.
    me_group = mp_rank(intra_group_comm)
    !
    ! ... search file
    rism1d_unit = find_free_unit()
    rism1d_file = TRIM(rism1d_file_base) // ".dat"
    exst = check_file_exst_1drism(TRIM(rism1d_file), ionode, ionode_id, intra_group_comm)
    !
    IF (.NOT. exst) THEN
      rism1d_file = TRIM(rism1d_file_base) // ".xml"
      exst = check_file_exst_1drism(TRIM(rism1d_file), ionode, ionode_id, intra_group_comm)
    END IF
    !
    IF (.NOT. exst) THEN
      CALL errore('read_1drism_xml', 'searching for ' // TRIM(rism1d_file), 10)
    END IF
    !
    ! ... open file
    IF (ionode) THEN
      CALL iotk_open_read(rism1d_unit, FILE=rism1d_file, IERR=ierr)
      CALL errore('read_1drism_xml', &
      & 'cannot open ' // TRIM(rism1d_file) // ' file for reading', ierr)
    END IF
    !
    ! ... read header
    IF (ionode) THEN
      CALL iotk_scan_begin(rism1d_unit, "_1D-RISM")
      !
      CALL iotk_scan_empty(rism1d_unit, "INFO", attr)
      CALL iotk_scan_attr(attr, "ngrid", ngrid_)
      CALL iotk_scan_attr(attr, "nsite", nsite_)
      !
      IF (ngrid /= ngrid_) THEN
        CALL errore('read_1drism_xml', 'number of grids do not match', 1)
      END IF
      !
      IF (nsite /= nsite_) THEN
        CALL errore('read_1drism_xml', 'number of sites do not match', 1)
      END IF
    END IF
    !
    ! ... find the index of the ionode
    io_group = 0
    IF (ionode) THEN
      io_group = me_group
    END IF
    CALL mp_sum(io_group, intra_group_comm)
    !
    ! ... read zvv for each site
    ALLOCATE(zvv1(ngrid))
    !
    DO isite = 1, nsite
      IF (ionode) THEN
        CALL iotk_scan_dat(rism1d_unit, "site" // iotk_index(isite), zvv1)
      END IF
      !
#if defined (__MPI)
      CALL MPI_SCATTERV(zvv1, npp, ipp, &
         & MPI_DOUBLE_PRECISION, zvv(1, isite), npp(me_group + 1), &
         & MPI_DOUBLE_PRECISION, io_group, intra_group_comm, ierr)
      !
      IF (ierr /= MPI_SUCCESS) THEN
        CALL errore('read_1drism_xml', 'error at MPI_SCATTERV', 1)
      END IF
#else
      zvv(1:ngrid, isite) = zvv1(1:ngrid)
#endif
    END DO
    !
    DEALLOCATE(zvv1)
    !
    ! ... close file
    IF (ionode) THEN
      CALL iotk_scan_end(rism1d_unit, "_1D-RISM")
      !
      CALL iotk_close_read(rism1d_unit)
    END IF
    !
  END SUBROUTINE read_1drism_xml
  !
  !------------------------------------------------------------------------
  FUNCTION check_file_exst_1drism(filename, ionode, ionode_id, intra_group_comm)
    !------------------------------------------------------------------------
    !
    !
    IMPLICIT NONE
    !
    LOGICAL :: check_file_exst_1drism
    !
    CHARACTER(LEN=*), INTENT(IN) :: filename
    LOGICAL,          INTENT(IN) :: ionode
    INTEGER,          INTENT(IN) :: ionode_id
    INTEGER,          INTENT(IN) :: intra_group_comm
    !
    LOGICAL :: lexists
    !
    IF (ionode) THEN
      !
      INQUIRE(FILE=TRIM(filename), EXIST=lexists)
      !
    END IF
    !
    CALL mp_bcast(lexists, ionode_id, intra_group_comm)
    !
    check_file_exst_1drism = lexists
    !
  END FUNCTION check_file_exst_1drism
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_3drism_xml(rism3d_file_base, zuv, name, &
                            & nsite, isite_start, isite_end, ecut, &
                            & nr1, nr2, nr3, nr1x, nr2x, ipp, npp, &
                            & ionode, intra_group_comm, inter_group_comm)
    !------------------------------------------------------------------------
    !
    ! ... Writs 3D-RISM's correlation function (R-space).
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: rism3d_file_base
    REAL(DP),         INTENT(IN) :: zuv(:,:)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER,          INTENT(IN) :: nsite
    INTEGER,          INTENT(IN) :: isite_start
    INTEGER,          INTENT(IN) :: isite_end
    REAL(DP),         INTENT(IN) :: ecut
    INTEGER,          INTENT(IN) :: nr1, nr2, nr3
    INTEGER,          INTENT(IN) :: nr1x, nr2x
    INTEGER,          INTENT(IN) :: ipp(:)
    INTEGER,          INTENT(IN) :: npp(:)
    LOGICAL,          INTENT(IN) :: ionode
    INTEGER,          INTENT(IN) :: inter_group_comm
    INTEGER,          INTENT(IN) :: intra_group_comm
    !
    INTEGER                 :: ip
    INTEGER                 :: i, j, k, kk
    INTEGER                 :: ierr
    INTEGER                 :: isite
    INTEGER                 :: iisite
    INTEGER                 :: rism3d_unit
    CHARACTER(LEN=256)      :: rism3d_file
    CHARACTER(LEN=10)       :: rism3d_extension
    CHARACTER(iotk_attlenx) :: attr
    INTEGER                 :: io_group
    INTEGER                 :: io_group_id
    INTEGER                 :: me_group
    INTEGER                 :: my_group_id
    INTEGER                 :: nproc_group
    INTEGER,  ALLOCATABLE   :: sowner(:)
    INTEGER,  ALLOCATABLE   :: kowner(:)
    REAL(DP), ALLOCATABLE   :: zuv_plane(:)
    !
    INTEGER, EXTERNAL       :: find_free_unit
    !
    ! ... allocate memory
    ALLOCATE(sowner(nsite))
    ALLOCATE(kowner(nr3))
    ALLOCATE(zuv_plane(nr1 * nr2))
    !
    ! ... get process info.
    me_group    = mp_rank(intra_group_comm)
    my_group_id = mp_rank(inter_group_comm)
    nproc_group = mp_size(intra_group_comm)
    !
    ! ... decide file name and unit
    rism3d_extension = '.dat'
    IF (.NOT. l3drism_binary) THEN
      rism3d_extension = '.xml'
    END IF
    !
    rism3d_file = TRIM(rism3d_file_base) // TRIM(rism3d_extension)
    rism3d_unit = find_free_unit()
    !
    ! ... open file
    IF (ionode) THEN
      CALL iotk_open_write(rism3d_unit, FILE=rism3d_file,  BINARY=l3drism_binary, IERR=ierr)
      CALL errore('write_rism3d_xml', &
      & 'cannot open ' // TRIM(rism3d_file) // ' file for writing', ierr)
    END IF
    !
    ! ... write header
    IF (ionode) THEN
      CALL iotk_write_begin(rism3d_unit, "_3D-RISM")
      !
      CALL iotk_write_attr(attr, "name", TRIM(name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "nsite", nsite)
      CALL iotk_write_attr(attr, "ecut", ecut)
      CALL iotk_write_attr(attr, "nr1", nr1)
      CALL iotk_write_attr(attr, "nr2", nr2)
      CALL iotk_write_attr(attr, "nr3", nr3)
      CALL iotk_write_empty(rism3d_unit, "INFO", attr)
    END IF
    !
    ! ... find the index of the group that will write zuv
    io_group_id = 0
    IF (ionode) THEN
      io_group_id = my_group_id
    END IF
    CALL mp_sum(io_group_id, intra_group_comm)
    CALL mp_sum(io_group_id, inter_group_comm)
    !
    ! ... find the index of the ionode within its own group
    io_group = 0
    IF (ionode) THEN
      io_group = me_group
    END IF
    CALL mp_sum(io_group, intra_group_comm)
    CALL mp_sum(io_group, inter_group_comm)
    !
    ! ... find out the owner of each solvent's site
    sowner = 0
    sowner(isite_start:isite_end) = my_group_id
    CALL mp_sum(sowner, inter_group_comm)
    !
    ! ... find out the owner of each "z" plane
    DO ip = 1, nproc_group
      kowner((ipp(ip) + 1):(ipp(ip) + npp(ip))) = ip - 1
    END DO
    !
    ! ... write zuv for each solvent's site
    DO isite = 1, nsite
      IF (ionode) THEN
        CALL iotk_write_begin(rism3d_unit, "site" // iotk_index(isite))
      END IF
      !
      IF (sowner(isite) == my_group_id) THEN
        iisite = isite - isite_start + 1
      ELSE
        iisite = -1
      END IF
      !
      ! ... write zuv for each "z" plane
      DO k = 1, nr3
        IF (sowner(isite) == my_group_id) THEN
          IF (kowner(k) == me_group) THEN
            kk = k - ipp(me_group + 1)
            DO j = 1, nr2
              DO i = 1, nr1
                zuv_plane(i + (j - 1) * nr1) = &
                & zuv(i + (j - 1) * nr1x + (kk - 1) * nr1x * nr2x, iisite)
              END DO
            END DO
          END IF
          !
          IF (kowner(k) /= io_group) THEN
            CALL mp_get(zuv_plane, zuv_plane, me_group, io_group, &
                      & kowner(k), k, intra_group_comm)
          END IF
        END IF
        !
        IF (sowner(isite) /= io_group_id) THEN
          CALL mp_get(zuv_plane, zuv_plane, my_group_id, io_group_id, &
                    & sowner(isite), isite, inter_group_comm)
        END IF
        !
        IF (ionode) THEN
          CALL iotk_write_dat(rism3d_unit, "z" // iotk_index(k), zuv_plane)
        END IF
      END DO
      !
      IF (ionode) THEN
        CALL iotk_write_end(rism3d_unit, "site" // iotk_index(isite))
      END IF
    END DO
    !
    ! ... close file
    IF (ionode) THEN
      CALL iotk_write_end(rism3d_unit, "_3D-RISM")
      !
      CALL iotk_close_write(rism3d_unit)
    END IF
    !
    ! ... deallocate memory
    DEALLOCATE(sowner)
    DEALLOCATE(kowner)
    DEALLOCATE(zuv_plane)
    !
  END SUBROUTINE write_3drism_xml
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_3drism_xml(rism3d_file_base, zuv, &
                           & nsite, isite_start, isite_end, ecut, &
                           & nr1, nr2, nr3, nr1x, nr2x, ipp, npp, &
                           & ionode, intra_group_comm, inter_group_comm)
    !------------------------------------------------------------------------
    !
    ! ... Reads 3D-RISM's correlation function (R-space).
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN)  :: rism3d_file_base
    REAL(DP),         INTENT(OUT) :: zuv(:,:)
    INTEGER,          INTENT(IN)  :: nsite
    INTEGER,          INTENT(IN)  :: isite_start
    INTEGER,          INTENT(IN)  :: isite_end
    REAL(DP),         INTENT(IN)  :: ecut
    INTEGER,          INTENT(IN)  :: nr1, nr2, nr3
    INTEGER,          INTENT(IN)  :: nr1x, nr2x
    INTEGER,          INTENT(IN)  :: ipp(:)
    INTEGER,          INTENT(IN)  :: npp(:)
    LOGICAL,          INTENT(IN)  :: ionode
    INTEGER,          INTENT(IN)  :: inter_group_comm
    INTEGER,          INTENT(IN)  :: intra_group_comm
    !
    INTEGER                 :: ip
    INTEGER                 :: i, j, k, kk
    INTEGER                 :: ierr
    INTEGER                 :: isite
    INTEGER                 :: iisite
    INTEGER                 :: rism3d_unit
    CHARACTER(LEN=256)      :: rism3d_file
    CHARACTER(iotk_attlenx) :: attr
    INTEGER                 :: nsite_
    REAL(DP)                :: ecut_
    INTEGER                 :: nr(3)
    INTEGER                 :: io_group
    INTEGER                 :: io_group_id
    INTEGER                 :: me_group
    INTEGER                 :: my_group_id
    INTEGER                 :: nproc_group
    LOGICAL                 :: exst
    INTEGER,  ALLOCATABLE   :: sowner(:)
    INTEGER,  ALLOCATABLE   :: kowner(:)
    REAL(DP), ALLOCATABLE   :: zuv_plane(:)
    !
    INTEGER, EXTERNAL       :: find_free_unit
    !
    ! ... allocate memory
    ALLOCATE(sowner(nsite))
    ALLOCATE(kowner(nr3))
    ALLOCATE(zuv_plane(nr1 * nr2))
    !
    ! ... get process info.
    me_group    = mp_rank(intra_group_comm)
    my_group_id = mp_rank(inter_group_comm)
    nproc_group = mp_size(intra_group_comm)
    !
    ! ... search file
    rism3d_unit = find_free_unit()
    rism3d_file = TRIM(rism3d_file_base) // ".dat"
    exst = check_file_exst(TRIM(rism3d_file))
    !
    IF (.NOT. exst) THEN
      rism3d_file = TRIM(rism3d_file_base) // ".xml"
      exst = check_file_exst(TRIM(rism3d_file))
    END IF
    !
    IF (.NOT. exst) THEN
      CALL errore('read_3drism_xml', 'searching for ' // TRIM(rism3d_file), 10)
    END IF
    !
    ! ... open file
    IF (ionode) THEN
      CALL iotk_open_read(rism3d_unit, FILE=rism3d_file, IERR=ierr)
      CALL errore('read_3drism_xml', &
      & 'cannot open ' // TRIM(rism3d_file) // ' file for reading', ierr)
    END IF
    !
    ! ... read header
    IF (ionode) THEN
      CALL iotk_scan_begin(rism3d_unit, "_3D-RISM")
      !
      CALL iotk_scan_empty(rism3d_unit, "INFO", attr)
      CALL iotk_scan_attr(attr, "nsite", nsite_)
      CALL iotk_scan_attr(attr, "ecut", ecut_)
      CALL iotk_scan_attr(attr, "nr1", nr(1))
      CALL iotk_scan_attr(attr, "nr2", nr(2))
      CALL iotk_scan_attr(attr, "nr3", nr(3))
      !
      IF (nsite /= nsite_) THEN
        CALL errore('read_3drism_xml', 'number of sites do not match', 1)
      END IF
      !
      IF (ABS(ecut - ecut_) > eps8) THEN
        CALL errore('read_3drism_xml', 'energy cutoff does not match', 1)
      END IF
      !
      IF (nr1 /= nr(1) .OR. nr2 /= nr(2) .OR. nr3 /= nr(3)) THEN
        CALL errore('read_3drism_xml', 'dimensions do not match', 1)
      END IF
    END IF
    !
    ! ... find the index of the group that will read zuv
    io_group_id = 0
    IF (ionode) THEN
      io_group_id = my_group_id
    END IF
    CALL mp_sum(io_group_id, intra_group_comm)
    CALL mp_sum(io_group_id, inter_group_comm)
    !
    ! ... find the index of the ionode within its own group
    io_group = 0
    IF (ionode) THEN
      io_group = me_group
    END IF
    CALL mp_sum(io_group, intra_group_comm)
    CALL mp_sum(io_group, inter_group_comm)
    !
    ! ... find out the owner of each solvent's site
    sowner = 0
    sowner(isite_start:isite_end) = my_group_id
    CALL mp_sum(sowner, inter_group_comm)
    !
    ! ... find out the owner of each "z" plane
    DO ip = 1, nproc_group
      kowner((ipp(ip) + 1):(ipp(ip) + npp(ip))) = ip - 1
    END DO
    !
    ! ... read zuv for each solvent's site
    DO isite = 1, nsite
      IF (ionode) THEN
        CALL iotk_scan_begin(rism3d_unit, "site" // iotk_index(isite))
      END IF
      !
      IF (sowner(isite) == my_group_id) THEN
        iisite = isite - isite_start + 1
      ELSE
        iisite = -1
      END IF
      !
      ! ... read zuv for each "z" plane
      DO k = 1, nr3
        IF (ionode) THEN
          CALL iotk_scan_dat(rism3d_unit, "z" // iotk_index(k), zuv_plane)
        END IF
        !
        IF (sowner(isite) /= io_group_id) THEN
          CALL mp_get(zuv_plane, zuv_plane, my_group_id, sowner(isite), &
                    & io_group_id, isite, inter_group_comm)
        END IF
        !
        IF (sowner(isite) == my_group_id) THEN
          IF (kowner(k) /= io_group) THEN
            CALL mp_get(zuv_plane, zuv_plane, me_group, kowner(k), &
                      & io_group, k, intra_group_comm)
          END IF
          !
          IF(kowner(k) == me_group) THEN
            kk = k - ipp(me_group + 1)
            DO j = 1, nr2
              DO i = 1, nr1
                zuv(i + (j - 1) * nr1x + (kk - 1 ) * nr1x * nr2x, iisite) &
                & = zuv_plane(i + (j - 1) * nr1)
              END DO
            END DO
          END IF
        END IF
        !
      END DO
      !
      IF (ionode) THEN
        CALL iotk_scan_end(rism3d_unit, "site" // iotk_index(isite))
      END IF
    END DO
    !
    ! ... close file
    IF (ionode) THEN
      CALL iotk_scan_end(rism3d_unit, "_3D-RISM")
      !
      CALL iotk_close_read(rism3d_unit)
    END IF
    !
    ! ... deallocate memory
    DEALLOCATE(sowner)
    DEALLOCATE(kowner)
    DEALLOCATE(zuv_plane)
    !
  END SUBROUTINE read_3drism_xml
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_lauerism_xml(rismlaue_file_base, zuv, name, &
                              & nsite, isite_start, isite_end, ecut, gamma_only, &
                              & lfft, ionode, intra_group_comm, inter_group_comm)
    !------------------------------------------------------------------------
    !
    ! ... Writs Laue-RISM's correlation function (Laue-rep., expanded Z-stick).
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),   INTENT(IN) :: rismlaue_file_base
    COMPLEX(DP),        INTENT(IN) :: zuv(:,:)
    CHARACTER(LEN=*),   INTENT(IN) :: name
    INTEGER,            INTENT(IN) :: nsite
    INTEGER,            INTENT(IN) :: isite_start
    INTEGER,            INTENT(IN) :: isite_end
    REAL(DP),           INTENT(IN) :: ecut
    LOGICAL,            INTENT(IN) :: gamma_only
    TYPE(lauefft_type), INTENT(IN) :: lfft
    LOGICAL,            INTENT(IN) :: ionode
    INTEGER,            INTENT(IN) :: inter_group_comm
    INTEGER,            INTENT(IN) :: intra_group_comm
    !
    INTEGER                  :: irz
    INTEGER                  :: igxy
    INTEGER                  :: jgxy1
    INTEGER                  :: jgxy2
    INTEGER                  :: igx, igy
    INTEGER                  :: mx, my
    INTEGER                  :: nr1, nr2, nr3
    INTEGER                  :: isign
    INTEGER                  :: ierr
    INTEGER                  :: isite
    INTEGER                  :: iisite
    REAL(DP)                 :: zreal
    REAL(DP)                 :: zimag
    INTEGER                  :: rismlaue_unit
    CHARACTER(LEN=256)       :: rismlaue_file
    CHARACTER(LEN=10)        :: rismlaue_extension
    CHARACTER(iotk_attlenx)  :: attr
    INTEGER                  :: io_group
    INTEGER                  :: io_group_id
    INTEGER                  :: me_group
    INTEGER                  :: my_group_id
    INTEGER,     ALLOCATABLE :: sowner(:)
    COMPLEX(DP), ALLOCATABLE :: zuv_site(:)
#if defined (__MPI)
    COMPLEX(DP), ALLOCATABLE :: zuv_tmp(:)
#endif
    !
    INTEGER, EXTERNAL        :: find_free_unit
    !
    ! ... set variables
    nr1 = lfft%dfft%nr1
    nr2 = lfft%dfft%nr2
    nr3 = lfft%nrz
    !
    ! ... allocate memory
    ALLOCATE(sowner(nsite))
    ALLOCATE(zuv_site(nr1 * nr2 * nr3))
    !
    ! ... get process info.
    me_group    = mp_rank(intra_group_comm)
    my_group_id = mp_rank(inter_group_comm)
    !
    ! ... decide file name and unit
    rismlaue_extension = '.dat'
    IF (.NOT. l3drism_binary) THEN
      rismlaue_extension = '.xml'
    END IF
    !
    rismlaue_file = TRIM(rismlaue_file_base) // TRIM(rismlaue_extension)
    rismlaue_unit = find_free_unit()
    !
    ! ... open file
    IF (ionode) THEN
      CALL iotk_open_write(rismlaue_unit, FILE=rismlaue_file,  BINARY=l3drism_binary, IERR=ierr)
      CALL errore('write_rismlaue_xml', &
      & 'cannot open ' // TRIM(rismlaue_file) // ' file for writing', ierr)
    END IF
    !
    ! ... write header
    IF (ionode) THEN
      CALL iotk_write_begin(rismlaue_unit, "LAUE-RISM")
      !
      CALL iotk_write_attr(attr, "name", TRIM(name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "nsite", nsite)
      CALL iotk_write_attr(attr, "ecut", ecut)
      CALL iotk_write_attr(attr, "nr1", nr1)
      CALL iotk_write_attr(attr, "nr2", nr2)
      CALL iotk_write_attr(attr, "nr3", nr3)
      CALL iotk_write_empty(rismlaue_unit, "INFO", attr)
    END IF
    !
    ! ... find the index of the group that will write zuv
    io_group_id = 0
    IF (ionode) THEN
      io_group_id = my_group_id
    END IF
    CALL mp_sum(io_group_id, intra_group_comm)
    CALL mp_sum(io_group_id, inter_group_comm)
    !
    ! ... find the index of the ionode within its own group
    io_group = 0
    IF (ionode) THEN
      io_group = me_group
    END IF
    CALL mp_sum(io_group, intra_group_comm)
    CALL mp_sum(io_group, inter_group_comm)
    !
    ! ... find out the owner of each solvent's site
    sowner = 0
    sowner(isite_start:isite_end) = my_group_id
    CALL mp_sum(sowner, inter_group_comm)
    !
    ! ... write zuv for each solvent's site
    DO isite = 1, nsite
      IF (sowner(isite) == my_group_id) THEN
        iisite = isite - isite_start + 1
      ELSE
        iisite = -1
      END IF
      !
      IF (sowner(isite) == my_group_id) THEN
        !
        CALL mp_barrier(intra_group_comm)
        !
        zuv_site = CMPLX(0.0_DP, 0.0_DP, kind=DP)
        !
        DO igxy = 1, lfft%ngxy
          isign = 1
       10 CONTINUE
          !
          mx  = isign * lfft%millxy(1, igxy)
          igx = mx + 1
          IF (igx < 1) THEN
            igx = igx + nr1
          END IF
          !
          my  = isign * lfft%millxy(2, igxy)
          igy = my + 1
          IF (igy < 1) THEN
            igy = igy + nr2
          END IF
          !
          jgxy1 = nr3 * (igxy - 1)
          jgxy2 = nr3 * (nr2 * (igx - 1) + (igy - 1))
          !
          DO irz = 1, nr3
            zreal = DBLE( zuv(irz + jgxy1, iisite))
            zimag = AIMAG(zuv(irz + jgxy1, iisite))
            zuv_site(irz + jgxy2) = CMPLX(zreal, DBLE(isign) * zimag, kind=DP)
          END DO
          !
          IF (gamma_only .AND. isign > 0) THEN
            isign = -1
            GOTO 10
          END IF
        END DO
        !
#if defined (__MPI)
#if defined (__FUJITSU)
        CALL mp_sum(zuv_site, intra_group_comm)
        !
#else
        IF (me_group == io_group) THEN
          ALLOCATE(zuv_tmp(nr1 * nr2 * nr3))
        ELSE
          ALLOCATE(zuv_tmp(1))
        END IF
        !
        CALL MPI_REDUCE(zuv_site(1), zuv_tmp(1), nr1 * nr2 * nr3, MPI_DOUBLE_COMPLEX, &
                      & MPI_SUM, io_group, intra_group_comm, ierr)
        !
        IF (ierr /= MPI_SUCCESS) THEN
          CALL errore('write_lauerism_xml', 'error at MPI_REDUCE', 1)
        END IF
        !
        IF (me_group == io_group) THEN
          zuv_site = zuv_tmp
        END IF
        DEALLOCATE(zuv_tmp)
        !
#endif
#endif
        !
      END IF
      !
      IF (sowner(isite) /= io_group_id) THEN
        IF (me_group == io_group) THEN
          !
          CALL mp_barrier(inter_group_comm)
          !
          CALL mp_get(zuv_site, zuv_site, my_group_id, io_group_id, &
                    & sowner(isite), isite, inter_group_comm)
        END IF
      END IF
      !
      IF (ionode) THEN
        CALL iotk_write_dat(rismlaue_unit, "site" // iotk_index(isite), zuv_site)
      END IF
    END DO
    !
    ! ... close file
    IF (ionode) THEN
      CALL iotk_write_end(rismlaue_unit, "LAUE-RISM")
      !
      CALL iotk_close_write(rismlaue_unit)
    END IF
    !
    ! ... deallocate memory
    DEALLOCATE(sowner)
    DEALLOCATE(zuv_site)
    !
  END SUBROUTINE write_lauerism_xml
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_lauerism_xml(rismlaue_file_base, zuv, &
                             & nsite, isite_start, isite_end, ecut, &
                             & lfft, ionode, intra_group_comm, inter_group_comm)
    !------------------------------------------------------------------------
    !
    ! ... Reads Laue-RISM's correlation function (Laue-rep., expanded Z-stick).
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),   INTENT(IN)  :: rismlaue_file_base
    COMPLEX(DP),        INTENT(OUT) :: zuv(:,:)
    INTEGER,            INTENT(IN)  :: nsite
    INTEGER,            INTENT(IN)  :: isite_start
    INTEGER,            INTENT(IN)  :: isite_end
    REAL(DP),           INTENT(IN)  :: ecut
    TYPE(lauefft_type), INTENT(IN)  :: lfft
    LOGICAL,            INTENT(IN)  :: ionode
    INTEGER,            INTENT(IN)  :: inter_group_comm
    INTEGER,            INTENT(IN)  :: intra_group_comm
    !
    INTEGER                  :: irz
    INTEGER                  :: igxy
    INTEGER                  :: jgxy1
    INTEGER                  :: jgxy2
    INTEGER                  :: igx, igy
    INTEGER                  :: mx, my
    INTEGER                  :: nr1, nr2, nr3
    INTEGER                  :: ierr
    INTEGER                  :: isite
    INTEGER                  :: iisite
    INTEGER                  :: rismlaue_unit
    CHARACTER(LEN=256)       :: rismlaue_file
    CHARACTER(iotk_attlenx)  :: attr
    INTEGER                  :: nsite_
    REAL(DP)                 :: ecut_
    INTEGER                  :: nr(3)
    INTEGER                  :: io_group
    INTEGER                  :: io_group_id
    INTEGER                  :: me_group
    INTEGER                  :: my_group_id
    LOGICAL                  :: exst
    INTEGER,     ALLOCATABLE :: sowner(:)
    COMPLEX(DP), ALLOCATABLE :: zuv_site(:)
    !
    INTEGER, EXTERNAL        :: find_free_unit
    !
    ! ... set variables
    nr1 = lfft%dfft%nr1
    nr2 = lfft%dfft%nr2
    nr3 = lfft%nrz
    !
    ! ... allocate memory
    ALLOCATE(sowner(nsite))
    ALLOCATE(zuv_site(nr1 * nr2 * nr3))
    !
    ! ... get process info.
    me_group    = mp_rank(intra_group_comm)
    my_group_id = mp_rank(inter_group_comm)
    !
    ! ... search file
    rismlaue_unit = find_free_unit()
    rismlaue_file = TRIM(rismlaue_file_base) // ".dat"
    exst = check_file_exst(TRIM(rismlaue_file))
    !
    IF (.NOT. exst) THEN
      rismlaue_file = TRIM(rismlaue_file_base) // ".xml"
      exst = check_file_exst(TRIM(rismlaue_file))
    END IF
    !
    IF (.NOT. exst) THEN
      CALL errore('read_lauerism_xml', 'searching for ' // TRIM(rismlaue_file), 10)
    END IF
    !
    ! ... open file
    IF (ionode) THEN
      CALL iotk_open_read(rismlaue_unit, FILE=rismlaue_file, IERR=ierr)
      CALL errore('read_lauerism_xml', &
      & 'cannot open ' // TRIM(rismlaue_file) // ' file for reading', ierr)
    END IF
    !
    ! ... read header
    IF (ionode) THEN
      CALL iotk_scan_begin(rismlaue_unit, "LAUE-RISM")
      !
      CALL iotk_scan_empty(rismlaue_unit, "INFO", attr)
      CALL iotk_scan_attr(attr, "nsite", nsite_)
      CALL iotk_scan_attr(attr, "ecut", ecut_)
      CALL iotk_scan_attr(attr, "nr1", nr(1))
      CALL iotk_scan_attr(attr, "nr2", nr(2))
      CALL iotk_scan_attr(attr, "nr3", nr(3))
      !
      IF (nsite /= nsite_) THEN
        CALL errore('read_lauerism_xml', 'number of sites do not match', 1)
      END IF
      !
      IF (ABS(ecut - ecut_) > eps8) THEN
        CALL errore('read_lauerism_xml', 'energy cutoff does not match', 1)
      END IF
      !
      IF (nr1 /= nr(1) .OR. nr2 /= nr(2) .OR. nr3 /= nr(3)) THEN
        CALL errore('read_lauerism_xml', 'dimensions do not match', 1)
      END IF
    END IF
    !
    ! ... find the index of the group that will read zuv
    io_group_id = 0
    IF (ionode) THEN
      io_group_id = my_group_id
    END IF
    CALL mp_sum(io_group_id, intra_group_comm)
    CALL mp_sum(io_group_id, inter_group_comm)
    !
    ! ... find the index of the ionode within its own group
    io_group = 0
    IF (ionode) THEN
      io_group = me_group
    END IF
    CALL mp_sum(io_group, intra_group_comm)
    CALL mp_sum(io_group, inter_group_comm)
    !
    ! ... find out the owner of each solvent's site
    sowner = 0
    sowner(isite_start:isite_end) = my_group_id
    CALL mp_sum(sowner, inter_group_comm)
    !
    ! ... read zuv for each solvent's site
    DO isite = 1, nsite
      IF (sowner(isite) == my_group_id) THEN
        iisite = isite - isite_start + 1
      ELSE
        iisite = -1
      END IF
      !
      IF (ionode) THEN
        CALL iotk_scan_dat(rismlaue_unit, "site" // iotk_index(isite), zuv_site)
      END IF
      !
      IF (my_group_id == io_group_id) THEN
        CALL mp_bcast(zuv_site, io_group, intra_group_comm)
      END IF
      !
      IF (sowner(isite) /= io_group_id) THEN
        CALL mp_get(zuv_site, zuv_site, my_group_id, sowner(isite), &
                  & io_group_id, isite, inter_group_comm)
      END IF
      !
      IF (sowner(isite) == my_group_id) THEN
        !
        DO igxy = 1, lfft%ngxy
          mx  = lfft%millxy(1, igxy)
          igx = mx + 1
          IF (igx < 1) THEN
            igx = igx + nr1
          END IF
          !
          my  = lfft%millxy(2, igxy)
          igy = my + 1
          IF (igy < 1) THEN
            igy = igy + nr2
          END IF
          !
          jgxy1 = nr3 * (igxy - 1)
          jgxy2 = nr3 * (nr2 * (igx - 1) + (igy - 1))
          !
          DO irz = 1, nr3
            zuv(irz + jgxy1, iisite) = zuv_site(irz + jgxy2)
          END DO
        END DO
        !
      END IF
    END DO
    !
    ! ... close file
    IF (ionode) THEN
      CALL iotk_scan_end(rismlaue_unit, "LAUE-RISM")
      !
      CALL iotk_close_read(rismlaue_unit)
    END IF
    !
    ! ... deallocate memory
    DEALLOCATE(sowner)
    DEALLOCATE(zuv_site)
    !
  END SUBROUTINE read_lauerism_xml
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_lauedipole_xml(rismlaue_file_base, zuv, name, &
                                & nsite, isite_start, isite_end, &
                                & ionode, intra_group_comm, inter_group_comm)
    !------------------------------------------------------------------------
    !
    ! ... Writs Laue-RISM's dipole function.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: rismlaue_file_base
    REAL(DP),         INTENT(IN) :: zuv(:)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER,          INTENT(IN) :: nsite
    INTEGER,          INTENT(IN) :: isite_start
    INTEGER,          INTENT(IN) :: isite_end
    LOGICAL,          INTENT(IN) :: ionode
    INTEGER,          INTENT(IN) :: inter_group_comm
    INTEGER,          INTENT(IN) :: intra_group_comm
    !
    INTEGER                  :: ierr
    INTEGER                  :: isite
    INTEGER                  :: iisite
    INTEGER                  :: rismlaue_unit
    CHARACTER(LEN=256)       :: rismlaue_file
    CHARACTER(LEN=10)        :: rismlaue_extension
    CHARACTER(iotk_attlenx)  :: attr
    INTEGER                  :: io_group
    INTEGER                  :: io_group_id
    INTEGER                  :: me_group
    INTEGER                  :: my_group_id
    INTEGER, ALLOCATABLE     :: sowner(:)
    REAL(DP)                 :: zuv_site
    REAL(DP)                 :: zuv_site_
    !
    INTEGER, EXTERNAL        :: find_free_unit
    !
    ! ... allocate memory
    ALLOCATE(sowner(nsite))
    !
    ! ... get process info.
    me_group    = mp_rank(intra_group_comm)
    my_group_id = mp_rank(inter_group_comm)
    !
    ! ... decide file name and unit
    !rismlaue_extension = '.dat'
    !IF (.NOT. l3drism_binary) THEN
    !  rismlaue_extension = '.xml'
    !END IF
    rismlaue_extension = '.xml'
    !
    rismlaue_file = TRIM(rismlaue_file_base) // TRIM(rismlaue_extension)
    rismlaue_unit = find_free_unit()
    !
    ! ... open file
    IF (ionode) THEN
      !CALL iotk_open_write(rismlaue_unit, FILE=rismlaue_file,  BINARY=l3drism_binary, IERR=ierr)
      CALL iotk_open_write(rismlaue_unit, FILE=rismlaue_file,  BINARY=.FALSE., IERR=ierr)
      CALL errore('write_lauedipole_xml', &
      & 'cannot open ' // TRIM(rismlaue_file) // ' file for writing', ierr)
    END IF
    !
    ! ... write header
    IF (ionode) THEN
      CALL iotk_write_begin(rismlaue_unit, "LAUE-RISM")
      !
      CALL iotk_write_attr(attr, "name", TRIM(name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "nsite", nsite)
      CALL iotk_write_empty(rismlaue_unit, "INFO", attr)
    END IF
    !
    ! ... find the index of the group that will write zuv
    io_group_id = 0
    IF (ionode) THEN
      io_group_id = my_group_id
    END IF
    CALL mp_sum(io_group_id, intra_group_comm)
    CALL mp_sum(io_group_id, inter_group_comm)
    !
    ! ... find the index of the ionode within its own group
    io_group = 0
    IF (ionode) THEN
      io_group = me_group
    END IF
    CALL mp_sum(io_group, intra_group_comm)
    CALL mp_sum(io_group, inter_group_comm)
    !
    ! ... find out the owner of each solvent's site
    sowner = 0
    sowner(isite_start:isite_end) = my_group_id
    CALL mp_sum(sowner, inter_group_comm)
    !
    ! ... write zuv for each solvent's site
    IF (me_group == io_group) THEN
      !
      DO isite = 1, nsite
        IF (sowner(isite) == my_group_id) THEN
          iisite = isite - isite_start + 1
        ELSE
          iisite = -1
        END IF
        !
        IF (sowner(isite) == my_group_id) THEN
          !
          zuv_site = zuv(iisite)
          !
        END IF
        !
        IF (sowner(isite) /= io_group_id) THEN
          !
          CALL mp_barrier(inter_group_comm)
          !
          CALL mp_get(zuv_site_, zuv_site, my_group_id, io_group_id, &
                    & sowner(isite), isite, inter_group_comm)
          !
          zuv_site = zuv_site_
          !
        END IF
        !
        IF (ionode) THEN
          CALL iotk_write_dat(rismlaue_unit, "site" // iotk_index(isite), zuv_site)
        END IF
      END DO
      !
    END IF
    !
    CALL mp_barrier(intra_group_comm)
    !
    ! ... close file
    IF (ionode) THEN
      CALL iotk_write_end(rismlaue_unit, "LAUE-RISM")
      !
      CALL iotk_close_write(rismlaue_unit)
    END IF
    !
    ! ... deallocate memory
    DEALLOCATE(sowner)
    !
  END SUBROUTINE write_lauedipole_xml
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_lauedipole_xml(rismlaue_file_base, zuv, &
                               & nsite, isite_start, isite_end, &
                               & ionode, intra_group_comm, inter_group_comm)
    !------------------------------------------------------------------------
    !
    ! ... Reads Laue-RISM's dipole function.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN)  :: rismlaue_file_base
    REAL(DP),         INTENT(OUT) :: zuv(:)
    INTEGER,          INTENT(IN)  :: nsite
    INTEGER,          INTENT(IN)  :: isite_start
    INTEGER,          INTENT(IN)  :: isite_end
    LOGICAL,          INTENT(IN)  :: ionode
    INTEGER,          INTENT(IN)  :: inter_group_comm
    INTEGER,          INTENT(IN)  :: intra_group_comm
    !
    INTEGER                 :: ierr
    INTEGER                 :: isite
    INTEGER                 :: iisite
    INTEGER                 :: rismlaue_unit
    CHARACTER(LEN=256)      :: rismlaue_file
    CHARACTER(iotk_attlenx) :: attr
    INTEGER                 :: nsite_
    INTEGER                 :: io_group
    INTEGER                 :: io_group_id
    INTEGER                 :: me_group
    INTEGER                 :: my_group_id
    LOGICAL                 :: exst
    INTEGER, ALLOCATABLE    :: sowner(:)
    REAL(DP)                :: zuv_site
    REAL(DP)                :: zuv_site_
    !
    INTEGER, EXTERNAL       :: find_free_unit
    !
    ! ... allocate memory
    ALLOCATE(sowner(nsite))
    !
    ! ... get process info.
    me_group    = mp_rank(intra_group_comm)
    my_group_id = mp_rank(inter_group_comm)
    !
    ! ... search file
    rismlaue_unit = find_free_unit()
    rismlaue_file = TRIM(rismlaue_file_base) // ".dat"
    exst = check_file_exst(TRIM(rismlaue_file))
    !
    IF (.NOT. exst) THEN
      rismlaue_file = TRIM(rismlaue_file_base) // ".xml"
      exst = check_file_exst(TRIM(rismlaue_file))
    END IF
    !
    IF (.NOT. exst) THEN
      CALL errore('read_lauedipole_xml', 'searching for ' // TRIM(rismlaue_file), 10)
    END IF
    !
    ! ... open file
    IF (ionode) THEN
      CALL iotk_open_read(rismlaue_unit, FILE=rismlaue_file, IERR=ierr)
      CALL errore('read_lauedipole_xml', &
      & 'cannot open ' // TRIM(rismlaue_file) // ' file for reading', ierr)
    END IF
    !
    ! ... read header
    IF (ionode) THEN
      CALL iotk_scan_begin(rismlaue_unit, "LAUE-RISM")
      !
      CALL iotk_scan_empty(rismlaue_unit, "INFO", attr)
      CALL iotk_scan_attr(attr, "nsite", nsite_)
      !
      IF (nsite /= nsite_) THEN
        CALL errore('read_lauedipole_xml', 'number of sites do not match', 1)
      END IF
    END IF
    !
    ! ... find the index of the group that will read zuv
    io_group_id = 0
    IF (ionode) THEN
      io_group_id = my_group_id
    END IF
    CALL mp_sum(io_group_id, intra_group_comm)
    CALL mp_sum(io_group_id, inter_group_comm)
    !
    ! ... find the index of the ionode within its own group
    io_group = 0
    IF (ionode) THEN
      io_group = me_group
    END IF
    CALL mp_sum(io_group, intra_group_comm)
    CALL mp_sum(io_group, inter_group_comm)
    !
    ! ... find out the owner of each solvent's site
    sowner = 0
    sowner(isite_start:isite_end) = my_group_id
    CALL mp_sum(sowner, inter_group_comm)
    !
    ! ... read zuv for each solvent's site
    IF (me_group == io_group) THEN
      !
      DO isite = 1, nsite
        IF (sowner(isite) == my_group_id) THEN
          iisite = isite - isite_start + 1
        ELSE
          iisite = -1
        END IF
        !
        IF (ionode) THEN
          CALL iotk_scan_dat(rismlaue_unit, "site" // iotk_index(isite), zuv_site)
        END IF
        !
        IF (sowner(isite) /= io_group_id) THEN
          !
          CALL mp_barrier(inter_group_comm)
          !
          CALL mp_get(zuv_site_, zuv_site, my_group_id, sowner(isite), &
                    & io_group_id, isite, inter_group_comm)
          !
          zuv_site = zuv_site_
          !
        END IF
        !
        IF (sowner(isite) == my_group_id) THEN
          !
          zuv(iisite) = zuv_site
          !
        END IF
      END DO
      !
    END IF
    !
    CALL mp_bcast(zuv, io_group, intra_group_comm)
    !
    ! ... close file
    IF (ionode) THEN
      CALL iotk_scan_end(rismlaue_unit, "LAUE-RISM")
      !
      CALL iotk_close_read(rismlaue_unit)
    END IF
    !
    ! ... deallocate memory
    DEALLOCATE(sowner)
    !
  END SUBROUTINE read_lauedipole_xml
  !
END MODULE xml_io_rism
