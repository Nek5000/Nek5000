PROGRAM cgnstonek

  USE CGNS
  USE module_global

  IMPLICIT NONE

  INTEGER(iprec) :: ier

  CHARACTER(132) :: error_message

  !
  CALL getarg(1,cgns_file)
  !
  IF (cgns_file == '') THEN
     !
     WRITE(6,'('' Input .cgns name: '')',ADVANCE = "NO")
     READ(5,'(a)') cgns_file
     !
  ENDIF
  cgns_len = INDEX (cgns_file, ' ') - 1
  cgns_file= cgns_file(1:cgns_len) //'.cgns'
  !
  ! open CGNS file for read
  !
  CALL cg_open_f(cgns_file, CG_MODE_READ, ifile, ier)

  IF (ier .NE. CG_OK) CALL cg_get_error_f(error_message)
  IF (ier .NE. CG_OK) CALL cg_error_exit_f
  !
  CALL cg_precision_f(ifile, pre, ier)
  CALL cg_version_f(ifile, version, ier)
  IF (ier .NE. 0) CALL cg_error_exit_f
  WRITE(6,'(''  cgns-version = '',f5.3)') version
  !
  CALL cg_nbases_f(ifile,nbases,ier)

  CALL read_cgns

  CALL write_re2_file
  !
  STOP
END PROGRAM cgnstonek

SUBROUTINE cgnstonek_exit (TEXT)
  !
  !********************************************************************
  !
  !     STOP PROGRAM
  !
  !
  USE module_global
  !
  IMPLICIT NONE
  !
  CHARACTER (*)  TEXT
  CHARACTER (132) :: file
  !
  !
  WRITE (6,*)
  WRITE (6,*)  TEXT
  WRITE (6,*)
  !
  !
  STOP

END SUBROUTINE cgnstonek_exit
