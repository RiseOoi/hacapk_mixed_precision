module m_HACApK_solve_methods
        use m_HACApK_base
        use m_HACApK_mine
        implicit none
        include 'mpif.h'
contains

! =============================================================================
! TRUE MIXED PRECISION DATA STRUCTURE H-MATRIX-VECTOR MULTIPLICATION!!!
! =============================================================================
subroutine  HACApK_adot_body_lfmtx_hyp_mixed(zau, st_leafmtxp_mixed, st_ctl, zu, nd)
  ! The main subroutine that computes the matrix-vector multiplication
  ! Input:
  !   st_leafmtxp: h-matrix
  !   st_ctl: a lot of parameters
  !   zu: source vector
  !   nd: dimension of source vector
  ! Returns:
  !   zau: result of matvecmul

  implicit none
  type(st_HACApK_leafmtxp_mixed) :: st_leafmtxp_mixed
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: zau(*)
  real*8 :: zu(*)

  real*8, dimension(:), allocatable :: zbut
  real*8, dimension(:), allocatable :: zaut

  integer*4, pointer :: lpmd(:), lnp(:), lsp(:), ltmp(:)
  1000 format(5(a, i10)/)
  2000 format(5(a, f10.4)/)

  integer, intent(in) :: nd

  ! double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
  integer :: mpinr, mpilog, nrank, icomm
  integer :: ip, i, k, t,  il, ill, it, itt
  integer :: nths, nthe, ith, ith1, ls, le
  integer :: ndl, ndt, kt, kt_64, kt_32, kt_tracker, ktmax, nlf, nstrtl, nstrtt
  integer :: calcOrder

  lpmd => st_ctl%lpmd(:)
  icomm = lpmd(1)
  nrank = lpmd(2)
  mpinr = lpmd(3)
  mpilog = lpmd(4)
  lnp(0:) => st_ctl%lnp
  lsp(0:) => st_ctl%lsp
  ltmp(0:) => st_ctl%lthr

  nlf = st_leafmtxp_mixed%nlf
  ktmax = st_leafmtxp_mixed%ktmax
  ith = omp_get_thread_num()  ! thread_num
  ith1 = ith+1
  nths = ltmp(ith)
  nthe = ltmp(ith1)-1
  allocate(zaut(nd))
  zaut(:) = 0.d0
  allocate(zbut(ktmax))
  ls = nd
  le = 1

  ! time_matvec_loop = 0.0d0

  do ip = nths, nthe
    ! We now have the coordinate of top left of submatrix (nstrtl & nstrtt)
    ! We can now use that information and row & column size of submatrix to
    ! loop over the whole thing and do something useful

    ndl = st_leafmtxp_mixed%st_lf(ip)%ndl
    ndt = st_leafmtxp_mixed%st_lf(ip)%ndt
    ! ns = ndl * ndt  ! useless
    nstrtl = st_leafmtxp_mixed%st_lf(ip)%nstrtl
    nstrtt = st_leafmtxp_mixed%st_lf(ip)%nstrtt

    if (nstrtl < ls) ls = nstrtl
    if (nstrtl + ndl - 1 > le) le = nstrtl + ndl - 1

    ! ltmtx == 1 means low-rank submatrix
    if (st_leafmtxp_mixed%st_lf(ip)%ltmtx == 1) then
      kt = st_leafmtxp_mixed%st_lf(ip)%kt
      kt_64 = st_leafmtxp_mixed%st_lf(ip)%kt_64
      kt_32 = st_leafmtxp_mixed%st_lf(ip)%kt_32
      zbut(1:kt) = 0.0

      ! a1 X source_vector part
      do il = 1, kt_64
        do it = 1, ndt
          itt = it + nstrtt - 1
          zbut(il) = zbut(il) + st_leafmtxp_mixed%st_lf(ip)%a1_64(it, il) * zu(itt)
        enddo
      enddo

      kt_tracker = 1
      do il = kt_64+1, kt
        do it = 1, ndt
          itt = it + nstrtt - 1
          zbut(il) = zbut(il) + st_leafmtxp_mixed%st_lf(ip)%a1_32(it, kt_tracker) * zu(itt)
        enddo
        kt_tracker = kt_tracker + 1
      enddo

      ! Diagonal multiplication part  ! TODO: WHETHER SHOULD I MERGE THIS DIAGONAL PART UP INTO THE TWO LOOPS ABOVE?
      do il = 1, kt
        zbut(il) = zbut(il) * st_leafmtxp_mixed%st_lf(ip)%diag(il)
      enddo

      ! a2 X halfway_vector part
      do il = 1, kt_64
        do it = 1, ndl
          ill = it + nstrtl - 1
          zaut(ill) = zaut(ill) + st_leafmtxp_mixed%st_lf(ip)%a2_64(it, il) * zbut(il)
        enddo
      enddo

      kt_tracker = 1
      do il = kt_64+1, kt
        do it = 1, ndl
          ill = it + nstrtl - 1
          zaut(ill) = zaut(ill) + st_leafmtxp_mixed%st_lf(ip)%a2_32(it, kt_tracker) * zbut(il)
        enddo
        kt_tracker = kt_tracker + 1
      enddo

    ! ltmtx == 2 means dense submatrix
    elseif (st_leafmtxp_mixed%st_lf(ip)%ltmtx == 2) then  ! (E1)
      do il = 1, ndl
        ill = il + nstrtl - 1
        do it = 1, ndt
          itt = it + nstrtt - 1
          zaut(ill) = zaut(ill) + st_leafmtxp_mixed%st_lf(ip)%a1_64(it, il) * zu(itt)
        enddo
      enddo
    endif

  enddo

  do il = ls, le
  !$omp atomic
    zau(il) = zau(il) + zaut(il)
  enddo

end subroutine HACApK_adot_body_lfmtx_hyp_mixed


! =============================================================================
! Scaled Leaves' Matvecmul
! =============================================================================
subroutine  HACApK_adot_body_lfmtx_hyp_scaled_64(zau, st_leafmtxp_scaled, st_ctl, zu, nd)
  ! The main subroutine that computes the matrix-vector multiplication
  ! Input:
  !   st_leafmtxp: h-matrix
  !   st_ctl: a lot of parameters
  !   zu: source vector
  !   nd: dimension of source vector
  ! Returns:
  !   zau: result of matvecmul

  implicit none
  type(st_HACApK_leafmtxp_scaled) :: st_leafmtxp_scaled
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: zau(*)
  real*8 :: zu(*)

  real*8, dimension(:), allocatable :: zbut
  real*8, dimension(:), allocatable :: zaut

  integer*4, pointer :: lpmd(:), lnp(:), lsp(:), ltmp(:)
  1000 format(5(a, i10)/)
  2000 format(5(a, f10.4)/)

  integer, intent(in) :: nd

  ! double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
  integer :: mpinr, mpilog, nrank, icomm
  integer :: ip, i, k, t,  il, ill, it, itt
  integer :: nths, nthe, ith, ith1, ls, le
  integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt
  integer :: calcOrder

  lpmd => st_ctl%lpmd(:)
  icomm = lpmd(1)
  nrank = lpmd(2)
  mpinr = lpmd(3)
  mpilog = lpmd(4)
  lnp(0:) => st_ctl%lnp
  lsp(0:) => st_ctl%lsp
  ltmp(0:) => st_ctl%lthr

  nlf = st_leafmtxp_scaled%nlf
  ktmax = st_leafmtxp_scaled%ktmax
  ith = omp_get_thread_num()  ! thread_num
  ith1 = ith+1
  nths = ltmp(ith)
  nthe = ltmp(ith1)-1
  allocate(zaut(nd))
  zaut(:) = 0.d0
  allocate(zbut(ktmax))
  ls = nd
  le = 1

  ! time_matvec_loop = 0.0d0

  do ip = nths, nthe
    ! We now have the coordinate of top left of submatrix (nstrtl & nstrtt)
    ! We can now use that information and row & column size of submatrix to
    ! loop over the whole thing and do something useful

    ndl = st_leafmtxp_scaled%st_lf(ip)%ndl
    ndt = st_leafmtxp_scaled%st_lf(ip)%ndt
    ! ns = ndl * ndt
    nstrtl = st_leafmtxp_scaled%st_lf(ip)%nstrtl
    nstrtt = st_leafmtxp_scaled%st_lf(ip)%nstrtt

    if (nstrtl < ls) ls = nstrtl
    if (nstrtl + ndl - 1 > le) le = nstrtl + ndl - 1

    ! time_matvec_loop_start = MPI_Wtime()

    ! ltmtx == 1 means low-rank submatrix
    if (st_leafmtxp_scaled%st_lf(ip)%ltmtx == 1) then
      ! Since we have decompose the submatrix into two matrices using
      ! Adaptive Cross-Approximation (ACA) techniques, instead of just
      ! size of (M X N), we now have two matrices of size (M X P) X (P X N)
      ! Here specifically, (ndl X ndt) => (ndl X kt) X (kt X ndt)
      kt = st_leafmtxp_scaled%st_lf(ip)%kt
      zbut(1:kt) = 0.0

      ! We now do a double-decomposed-matrices-vector multiplication
      ! (ndl X kt) X (kt X ndt) X (ndl)

      ! The loop below is doing a (kt X ndt) X (ndl X 1) matvec_muliplication
      ! The result is a (kt X 1) vector
      do il = 1, kt
        do it = 1 , ndt
          itt = it + nstrtt - 1
          zbut(il) = zbut(il) + st_leafmtxp_scaled%st_lf(ip)%a1(it,il) * zu(itt)
        enddo
        ! Diagonal multiplication part
        zbut(il) = zbut(il) * st_leafmtxp_scaled%st_lf(ip)%diag(il)
      enddo

      ! outer diagonal test
      ! do il = 1, kt
      !   zbut(il) = zbut(il) * st_leafmtxp_scaled%st_lf(ip)%diag(il)
      ! enddo

      ! The loop below is doing a (ndl X kt) X (kt X 1) matvec_muliplication
      ! The result is a (ndl X 1) vector
      do il = 1, kt
        do it = 1, ndl
          ill = it + nstrtl - 1
          zaut(ill) = zaut(ill) + st_leafmtxp_scaled%st_lf(ip)%a2(it,il) * zbut(il)
        enddo
      enddo

    ! ltmtx == 2 means dense submatrix
    elseif (st_leafmtxp_scaled%st_lf(ip)%ltmtx == 2) then  ! (E1)
      do il = 1, ndl
        ill = il + nstrtl - 1
        do it = 1, ndt
          itt = it + nstrtt - 1
          zaut(ill) = zaut(ill) + st_leafmtxp_scaled%st_lf(ip)%a1(it,il) * zu(itt)
        enddo
      enddo
    endif

  ! time_matvec_loop_end = MPI_Wtime()
  ! time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start

  enddo

  ! print *, "Time of matvec loop:", time_matvec_loop

  do il = ls, le
  !$omp atomic
    zau(il) = zau(il) + zaut(il)
  enddo

end subroutine HACApK_adot_body_lfmtx_hyp_scaled_64


subroutine  HACApK_adot_body_lfmtx_hyp_scaled_32(zau, st_leafmtxp_scaled_32, st_ctl, zu_32, nd)
  ! The main subroutine that computes the matrix-vector multiplication
  ! Input:
  !   st_leafmtxp: h-matrix
  !   st_ctl: a lot of parameters
  !   zu_32: source vector
  !   nd: dimension of source vector
  ! Returns:
  !   zau: result of matvecmul

  implicit none
  type(st_HACApK_leafmtxp_scaled_32) :: st_leafmtxp_scaled_32
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: zau(*)
  real*4 :: zu_32(*)

  real*8, dimension(:), allocatable :: zbut
  real*8, dimension(:), allocatable :: zaut

  integer*4, pointer :: lpmd(:), lnp(:), lsp(:), ltmp(:)
  1000 format(5(a, i10)/)
  2000 format(5(a, f10.4)/)

  integer, intent(in) :: nd

  ! double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
  integer :: mpinr, mpilog, nrank, icomm
  integer :: ip, i, k, t,  il, ill, it, itt
  integer :: nths, nthe, ith, ith1, ls, le
  integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt
  integer :: calcOrder

  lpmd => st_ctl%lpmd(:)
  icomm = lpmd(1)
  nrank = lpmd(2)
  mpinr = lpmd(3)
  mpilog = lpmd(4)
  lnp(0:) => st_ctl%lnp
  lsp(0:) => st_ctl%lsp
  ltmp(0:) => st_ctl%lthr

  nlf = st_leafmtxp_scaled_32%nlf
  ktmax = st_leafmtxp_scaled_32%ktmax
  ith = omp_get_thread_num()  ! thread_num
  ith1 = ith+1
  nths = ltmp(ith)
  nthe = ltmp(ith1)-1
  allocate(zaut(nd))
  zaut(:) = 0.d0
  allocate(zbut(ktmax))
  ls = nd
  le = 1

  ! time_matvec_loop = 0.0d0

  do ip = nths, nthe
    ! We now have the coordinate of top left of submatrix (nstrtl & nstrtt)
    ! We can now use that information and row & column size of submatrix to
    ! loop over the whole thing and do something useful

    ndl = st_leafmtxp_scaled_32%st_lf(ip)%ndl
    ndt = st_leafmtxp_scaled_32%st_lf(ip)%ndt
    ! ns = ndl * ndt
    nstrtl = st_leafmtxp_scaled_32%st_lf(ip)%nstrtl
    nstrtt = st_leafmtxp_scaled_32%st_lf(ip)%nstrtt

    if (nstrtl < ls) ls = nstrtl
    if (nstrtl + ndl - 1 > le) le = nstrtl + ndl - 1

    ! time_matvec_loop_start = MPI_Wtime()

    ! ltmtx == 1 means low-rank submatrix
    if (st_leafmtxp_scaled_32%st_lf(ip)%ltmtx == 1) then
      ! Since we have decompose the submatrix into two matrices using
      ! Adaptive Cross-Approximation (ACA) techniques, instead of just
      ! size of (M X N), we now have two matrices of size (M X P) X (P X N)
      ! Here specifically, (ndl X ndt) => (ndl X kt) X (kt X ndt)
      kt = st_leafmtxp_scaled_32%st_lf(ip)%kt
      zbut(1:kt) = 0.0

      ! We now do a double-decomposed-matrices-vector multiplication
      ! (ndl X kt) X (kt X ndt) X (ndl)

      ! The loop below is doing a (kt X ndt) X (ndl X 1) matvec_muliplication
      ! The result is a (kt X 1) vector
      do il = 1, kt
        do it = 1 , ndt
          itt = it + nstrtt - 1
          zbut(il) = zbut(il) + st_leafmtxp_scaled_32%st_lf(ip)%a1(it,il) * zu_32(itt)
        enddo
        ! Diagonal multiplication part
        zbut(il) = zbut(il) * st_leafmtxp_scaled_32%st_lf(ip)%diag(il)
      enddo

      ! outer diag test
      ! do il = 1, kt
      !   zbut(il) = zbut(il) * st_leafmtxp_scaled_32%st_lf(ip)%diag(il)
      ! enddo

      ! The loop below is doing a (ndl X kt) X (kt X 1) matvec_muliplication
      ! The result is a (ndl X 1) vector
      do il = 1, kt
        do it = 1, ndl
          ill = it + nstrtl - 1
          zaut(ill) = zaut(ill) + st_leafmtxp_scaled_32%st_lf(ip)%a2(it,il) * zbut(il)
        enddo
      enddo

    ! ltmtx == 2 means dense submatrix
    elseif (st_leafmtxp_scaled_32%st_lf(ip)%ltmtx == 2) then  ! (E1)
      do il = 1, ndl
        ill = il + nstrtl - 1
        do it = 1, ndt
          itt = it + nstrtt - 1
          zaut(ill) = zaut(ill) + st_leafmtxp_scaled_32%st_lf(ip)%a1(it,il) * zu_32(itt)
        enddo
      enddo
    endif

  ! time_matvec_loop_end = MPI_Wtime()
  ! time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start

  enddo

  ! print *, "Time of matvec loop:", time_matvec_loop

  do il = ls, le
  !$omp atomic
    zau(il) = zau(il) + zaut(il)
  enddo

end subroutine HACApK_adot_body_lfmtx_hyp_scaled_32


subroutine  HACApK_adot_body_lfmtx_hyp_scaled_mixed(zau, st_leafmtxp_scaled_32, st_ctl, zu, nd)
  ! The main subroutine that computes the matrix-vector multiplication
  ! Input:
  !   st_leafmtxp: h-matrix
  !   st_ctl: a lot of parameters
  !   zu: source vector
  !   nd: dimension of source vector
  ! Returns:
  !   zau: result of matvecmul

  implicit none
  type(st_HACApK_leafmtxp_scaled_32) :: st_leafmtxp_scaled_32
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: zau(*)
  real*8 :: zu(*)

  real*8, dimension(:), allocatable :: zbut
  real*8, dimension(:), allocatable :: zaut

  integer*4, pointer :: lpmd(:), lnp(:), lsp(:), ltmp(:)
  1000 format(5(a, i10)/)
  2000 format(5(a, f10.4)/)

  integer, intent(in) :: nd

  ! double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
  integer :: mpinr, mpilog, nrank, icomm
  integer :: ip, i, k, t,  il, ill, it, itt
  integer :: nths, nthe, ith, ith1, ls, le
  integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt
  integer :: calcOrder

  lpmd => st_ctl%lpmd(:)
  icomm = lpmd(1)
  nrank = lpmd(2)
  mpinr = lpmd(3)
  mpilog = lpmd(4)
  lnp(0:) => st_ctl%lnp
  lsp(0:) => st_ctl%lsp
  ltmp(0:) => st_ctl%lthr

  nlf = st_leafmtxp_scaled_32%nlf
  ktmax = st_leafmtxp_scaled_32%ktmax
  ith = omp_get_thread_num()  ! thread_num
  ith1 = ith+1
  nths = ltmp(ith)
  nthe = ltmp(ith1)-1
  allocate(zaut(nd))
  zaut(:) = 0.d0
  allocate(zbut(ktmax))
  ls = nd
  le = 1

  ! time_matvec_loop = 0.0d0

  do ip = nths, nthe
    ! We now have the coordinate of top left of submatrix (nstrtl & nstrtt)
    ! We can now use that information and row & column size of submatrix to
    ! loop over the whole thing and do something useful

    ndl = st_leafmtxp_scaled_32%st_lf(ip)%ndl
    ndt = st_leafmtxp_scaled_32%st_lf(ip)%ndt
    ! ns = ndl * ndt
    nstrtl = st_leafmtxp_scaled_32%st_lf(ip)%nstrtl
    nstrtt = st_leafmtxp_scaled_32%st_lf(ip)%nstrtt

    if (nstrtl < ls) ls = nstrtl
    if (nstrtl + ndl - 1 > le) le = nstrtl + ndl - 1

    ! time_matvec_loop_start = MPI_Wtime()

    ! ltmtx == 1 means low-rank submatrix
    if (st_leafmtxp_scaled_32%st_lf(ip)%ltmtx == 1) then
      ! Since we have decompose the submatrix into two matrices using
      ! Adaptive Cross-Approximation (ACA) techniques, instead of just
      ! size of (M X N), we now have two matrices of size (M X P) X (P X N)
      ! Here specifically, (ndl X ndt) => (ndl X kt) X (kt X ndt)
      kt = st_leafmtxp_scaled_32%st_lf(ip)%kt
      zbut(1:kt) = 0.0

      ! We now do a double-decomposed-matrices-vector multiplication
      ! (ndl X kt) X (kt X ndt) X (ndl)

      ! The loop below is doing a (kt X ndt) X (ndl X 1) matvec_muliplication
      ! The result is a (kt X 1) vector
      do il = 1, kt
        do it = 1 , ndt
          itt = it + nstrtt - 1
          zbut(il) = zbut(il) + st_leafmtxp_scaled_32%st_lf(ip)%a1(it,il) * zu(itt)
        enddo
        ! Diagonal multiplication part
        zbut(il) = zbut(il) * st_leafmtxp_scaled_32%st_lf(ip)%diag(il)
      enddo

      ! outer diagonal test
      ! do il = 1, kt
      !   zbut(il) = zbut(il) * st_leafmtxp_scaled_32%st_lf(ip)%diag(il)
      ! enddo

      ! The loop below is doing a (ndl X kt) X (kt X 1) matvec_muliplication
      ! The result is a (ndl X 1) vector
      do il = 1, kt
        do it = 1, ndl
          ill = it + nstrtl - 1
          zaut(ill) = zaut(ill) + st_leafmtxp_scaled_32%st_lf(ip)%a2(it,il) * zbut(il)
        enddo
      enddo

    ! ltmtx == 2 means dense submatrix
    elseif (st_leafmtxp_scaled_32%st_lf(ip)%ltmtx == 2) then  ! (E1)
      do il = 1, ndl
        ill = il + nstrtl - 1
        do it = 1, ndt
          itt = it + nstrtt - 1
          zaut(ill) = zaut(ill) + st_leafmtxp_scaled_32%st_lf(ip)%a1(it,il) * zu(itt)
        enddo
      enddo
    endif

  ! time_matvec_loop_end = MPI_Wtime()
  ! time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start

  enddo

  ! print *, "Time of matvec loop:", time_matvec_loop

  do il = ls, le
  !$omp atomic
    zau(il) = zau(il) + zaut(il)
  enddo

end subroutine HACApK_adot_body_lfmtx_hyp_scaled_mixed
! =============================================================================
! Scaling methods stop here
! =============================================================================


! =============================================================================
! Rise 32 x 64
! =============================================================================
subroutine  HACApK_adot_body_lfmtx_hyp_rise_32_64(zau, st_leafmtxp_32, st_ctl, zu, nd)
  ! The main subroutine that computes the matrix-vector multiplication
  ! Input:
  !   st_leafmtxp: h-matrix
  !   st_ctl: a lot of parameters
  !   zu: source vector
  !   nd: dimension of source vector
  ! Returns:
  !   zau: result of matvecmul

  implicit none
  type(st_HACApK_leafmtxp_32) :: st_leafmtxp_32
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: zau(*)
  real*8 :: zu(*)  ! assume zu as 64 bit for now

  real*4, dimension(:), allocatable :: zbut_32
  real*4, dimension(:), allocatable :: zaut_32

  integer*4, pointer :: lpmd(:), lnp(:), lsp(:), ltmp(:)
  1000 format(5(a, i10)/)
  2000 format(5(a, f10.4)/)

  integer, intent(in) :: nd

  ! double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
  integer :: mpinr, mpilog, nrank, icomm
  integer :: ip, i, k, t,  il, ill, it, itt
  integer :: nths, nthe, ith, ith1, ls, le
  integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt
  integer :: calcOrder

  lpmd => st_ctl%lpmd(:)
  icomm = lpmd(1)
  nrank = lpmd(2)
  mpinr = lpmd(3)
  mpilog = lpmd(4)
  lnp(0:) => st_ctl%lnp
  lsp(0:) => st_ctl%lsp
  ltmp(0:) => st_ctl%lthr

  nlf = st_leafmtxp_32%nlf
  ktmax = st_leafmtxp_32%ktmax
  ith = omp_get_thread_num()  ! thread_num
  ith1 = ith+1
  nths = ltmp(ith)
  nthe = ltmp(ith1)-1
  allocate(zaut_32(nd))
  zaut_32(:) = 0.d0
  allocate(zbut_32(ktmax))
  ls = nd
  le = 1

  ! time_matvec_loop = 0.0d0

  do ip = nths, nthe
    ! We now have the coordinate of top left of submatrix (nstrtl & nstrtt)
    ! We can now use that information and row & column size of submatrix to
    ! loop over the whole thing and do something useful
    ndl = st_leafmtxp_32%st_lf(ip)%ndl
    ndt = st_leafmtxp_32%st_lf(ip)%ndt
    ! ns = ndl * ndt  ! (E1)
    nstrtl = st_leafmtxp_32%st_lf(ip)%nstrtl
    nstrtt = st_leafmtxp_32%st_lf(ip)%nstrtt

    if (nstrtl < ls) ls = nstrtl
    if (nstrtl + ndl - 1 > le) le = nstrtl + ndl - 1

    ! time_matvec_loop_start = MPI_Wtime()

    ! ltmtx == 1 means low-rank submatrix
    if (st_leafmtxp_32%st_lf(ip)%ltmtx == 1) then
      ! Since we have decompose the submatrix into two matrices using
      ! Adaptive Cross-Approximation (ACA) techniques, instead of just
      ! size of (M X N), we now have two matrices of size (M X P) X (P X N)
      ! Here specifically, (ndl X ndt) => (ndl X kt) X (kt X ndt)
      kt = st_leafmtxp_32%st_lf(ip)%kt
      zbut_32(1:kt) = 0.0

      ! We now do a double-decomposed-matrices-vector multiplication
      ! (ndl X kt) X (kt X ndt) X (ndl)

      ! The loop below is doing a (kt X ndt) X (ndl X 1) matvec_muliplication
      ! The result is a (kt X 1) vector
      do il = 1, kt
        do it = 1 , ndt
          itt = it + nstrtt - 1
          zbut_32(il) = zbut_32(il) + st_leafmtxp_32%st_lf(ip)%a1(it,il) * zu(itt)
        enddo
      enddo

      ! The loop below is doing a (ndl X kt) X (kt X 1) matvec_muliplication
      ! The result is a (ndl X 1) vector
      do il = 1, kt
        do it = 1, ndl
          ill = it + nstrtl - 1
          zaut_32(ill) = zaut_32(ill) + st_leafmtxp_32%st_lf(ip)%a2(it,il) * zbut_32(il)
        enddo
      enddo

    ! ltmtx == 2 means dense submatrix
    elseif (st_leafmtxp_32%st_lf(ip)%ltmtx == 2) then  ! (E1)
      do il = 1, ndl
        ill = il + nstrtl - 1
        do it = 1, ndt
          itt = it + nstrtt - 1
          zaut_32(ill) = zaut_32(ill) + st_leafmtxp_32%st_lf(ip)%a1(it,il) * zu(itt)
        enddo
      enddo
    endif

  ! time_matvec_loop_end = MPI_Wtime()
  ! time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start

  enddo

  ! print *, "Time of matvec loop:", time_matvec_loop

  do il = ls, le
  !$omp atomic
    zau(il) = zau(il) + zaut_32(il)
  enddo

end subroutine HACApK_adot_body_lfmtx_hyp_rise_32_64
! =============================================================================
! Rise 23 method stops here
! =============================================================================


! =============================================================================
! Rise 24  ! (E1)
! The goal of this method is to allow lower precision submatrices
! and check if the speed of matrix-vector computation increases.
! (32-bit declaration outside this method)
! =============================================================================
subroutine  HACApK_adot_body_lfmtx_hyp_rise_32(zau, st_leafmtxp_32, st_ctl, zu_32, nd)
  ! The main subroutine that computes the matrix-vector multiplication
  ! Input:
  !   st_leafmtxp: h-matrix
  !   st_ctl: a lot of parameters
  !   zu: source vector
  !   nd: dimension of source vector
  ! Returns:
  !   zau: result of matvecmul

  implicit none
  type(st_HACApK_leafmtxp_32) :: st_leafmtxp_32
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: zau(*)
  real*4 :: zu_32(*)
  ! real*8 :: zu_32(*)

  real*4, dimension(:), allocatable :: zbut_32
  real*4, dimension(:), allocatable :: zaut_32

  integer*4, pointer :: lpmd(:), lnp(:), lsp(:), ltmp(:)
  1000 format(5(a, i10)/)
  2000 format(5(a, f10.4)/)

  integer, intent(in) :: nd

  ! double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
  integer :: mpinr, mpilog, nrank, icomm
  integer :: ip, i, k, t,  il, ill, it, itt
  integer :: nths, nthe, ith, ith1, ls, le
  integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt
  integer :: calcOrder

  lpmd => st_ctl%lpmd(:)
  icomm = lpmd(1)
  nrank = lpmd(2)
  mpinr = lpmd(3)
  mpilog = lpmd(4)
  lnp(0:) => st_ctl%lnp
  lsp(0:) => st_ctl%lsp
  ltmp(0:) => st_ctl%lthr

  nlf = st_leafmtxp_32%nlf
  ktmax = st_leafmtxp_32%ktmax
  ith = omp_get_thread_num()  ! thread_num
  ith1 = ith+1
  nths = ltmp(ith)
  nthe = ltmp(ith1)-1
  allocate(zaut_32(nd))
  zaut_32(:) = 0.d0
  allocate(zbut_32(ktmax))
  ls = nd
  le = 1

  ! time_matvec_loop = 0.0d0

  do ip = nths, nthe
    ! We now have the coordinate of top left of submatrix (nstrtl & nstrtt)
    ! We can now use that information and row & column size of submatrix to
    ! loop over the whole thing and do something useful
    ndl = st_leafmtxp_32%st_lf(ip)%ndl
    ndt = st_leafmtxp_32%st_lf(ip)%ndt
    ! ns = ndl * ndt  ! (E1)
    nstrtl = st_leafmtxp_32%st_lf(ip)%nstrtl
    nstrtt = st_leafmtxp_32%st_lf(ip)%nstrtt

    if (nstrtl < ls) ls = nstrtl
    if (nstrtl + ndl - 1 > le) le = nstrtl + ndl - 1

    ! time_matvec_loop_start = MPI_Wtime()

    ! ltmtx == 1 means low-rank submatrix
    if (st_leafmtxp_32%st_lf(ip)%ltmtx == 1) then
      ! Since we have decompose the submatrix into two matrices using
      ! Adaptive Cross-Approximation (ACA) techniques, instead of just
      ! size of (M X N), we now have two matrices of size (M X P) X (P X N)
      ! Here specifically, (ndl X ndt) => (ndl X kt) X (kt X ndt)
      kt = st_leafmtxp_32%st_lf(ip)%kt
      zbut_32(1:kt) = 0.0

      ! We now do a double-decomposed-matrices-vector multiplication
      ! (ndl X kt) X (kt X ndt) X (ndl)

      ! The loop below is doing a (kt X ndt) X (ndl X 1) matvec_muliplication
      ! The result is a (kt X 1) vector
      do il = 1, kt
        do it = 1 , ndt
          itt = it + nstrtt - 1
          zbut_32(il) = zbut_32(il) + st_leafmtxp_32%st_lf(ip)%a1(it,il) * zu_32(itt)
        enddo
      enddo

      ! The loop below is doing a (ndl X kt) X (kt X 1) matvec_muliplication
      ! The result is a (ndl X 1) vector
      do il = 1, kt
        do it = 1, ndl
          ill = it + nstrtl - 1
          zaut_32(ill) = zaut_32(ill) + st_leafmtxp_32%st_lf(ip)%a2(it,il) * zbut_32(il)
        enddo
      enddo

    ! ltmtx == 2 means dense submatrix
    elseif (st_leafmtxp_32%st_lf(ip)%ltmtx == 2) then  ! (E1)
      do il = 1, ndl
        ill = il + nstrtl - 1
        do it = 1, ndt
          itt = it + nstrtt - 1
          zaut_32(ill) = zaut_32(ill) + st_leafmtxp_32%st_lf(ip)%a1(it,il) * zu_32(itt)
        enddo
      enddo
    endif

  ! time_matvec_loop_end = MPI_Wtime()
  ! time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start

  enddo

  ! print *, "Time of matvec loop:", time_matvec_loop

  do il = ls, le
  !$omp atomic
    zau(il) = zau(il) + zaut_32(il)
  enddo

end subroutine HACApK_adot_body_lfmtx_hyp_rise_32
! =============================================================================
! Rise 24 method stops here
! =============================================================================


! =============================================================================
! Rise 25  ! (E1)
! The goal of this method is to allow lower precision submatrices
! and check if the speed of matrix-vector computation increases.
! (32-bit declaration outside this method)
! =============================================================================
subroutine  HACApK_adot_body_lfmtx_hyp_rise_64(zau, st_leafmtxp, st_ctl, zu, nd)
  ! The main subroutine that computes the matrix-vector multiplication
  ! Input:
  !   st_leafmtxp: h-matrix
  !   st_ctl: a lot of parameters
  !   zu: source vector
  !   nd: dimension of source vector
  ! Returns:
  !   zau: result of matvecmul

  implicit none
  type(st_HACApK_leafmtxp) :: st_leafmtxp
  type(st_HACApK_lcontrol) :: st_ctl
  real*8 :: zau(*)
  real*8 :: zu(*)

  real*8, dimension(:), allocatable :: zbut
  real*8, dimension(:), allocatable :: zaut

  integer*4, pointer :: lpmd(:), lnp(:), lsp(:), ltmp(:)
  1000 format(5(a, i10)/)
  2000 format(5(a, f10.4)/)

  integer, intent(in) :: nd

  ! double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
  integer :: mpinr, mpilog, nrank, icomm
  integer :: ip, i, k, t,  il, ill, it, itt
  integer :: nths, nthe, ith, ith1, ls, le
  integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt
  integer :: calcOrder

  lpmd => st_ctl%lpmd(:)
  icomm = lpmd(1)
  nrank = lpmd(2)
  mpinr = lpmd(3)
  mpilog = lpmd(4)
  lnp(0:) => st_ctl%lnp
  lsp(0:) => st_ctl%lsp
  ltmp(0:) => st_ctl%lthr

  nlf = st_leafmtxp%nlf
  ktmax = st_leafmtxp%ktmax
  ith = omp_get_thread_num()  ! thread_num
  ith1 = ith+1
  nths = ltmp(ith)
  nthe = ltmp(ith1)-1
  allocate(zaut(nd))
  zaut(:) = 0.d0
  allocate(zbut(ktmax))
  ls = nd
  le = 1

  ! time_matvec_loop = 0.0d0

  do ip = nths, nthe
    ! We now have the coordinate of top left of submatrix (nstrtl & nstrtt)
    ! We can now use that information and row & column size of submatrix to
    ! loop over the whole thing and do something useful

    ndl = st_leafmtxp%st_lf(ip)%ndl
    ndt = st_leafmtxp%st_lf(ip)%ndt
    ! ns = ndl * ndt
    nstrtl = st_leafmtxp%st_lf(ip)%nstrtl
    nstrtt = st_leafmtxp%st_lf(ip)%nstrtt

    if (nstrtl < ls) ls = nstrtl
    if (nstrtl + ndl - 1 > le) le = nstrtl + ndl - 1

    ! time_matvec_loop_start = MPI_Wtime()

    ! ltmtx == 1 means low-rank submatrix
    if (st_leafmtxp%st_lf(ip)%ltmtx == 1) then
      ! Since we have decompose the submatrix into two matrices using
      ! Adaptive Cross-Approximation (ACA) techniques, instead of just
      ! size of (M X N), we now have two matrices of size (M X P) X (P X N)
      ! Here specifically, (ndl X ndt) => (ndl X kt) X (kt X ndt)
      kt = st_leafmtxp%st_lf(ip)%kt
      zbut(1:kt) = 0.0

      ! We now do a double-decomposed-matrices-vector multiplication
      ! (ndl X kt) X (kt X ndt) X (ndl)

      ! The loop below is doing a (kt X ndt) X (ndl X 1) matvec_muliplication
      ! The result is a (kt X 1) vector
      do il = 1, kt
        do it = 1 , ndt
          itt = it + nstrtt - 1
          zbut(il) = zbut(il) + st_leafmtxp%st_lf(ip)%a1(it,il) * zu(itt)
        enddo
      enddo

      ! The loop below is doing a (ndl X kt) X (kt X 1) matvec_muliplication
      ! The result is a (ndl X 1) vector
      do il = 1, kt
        do it = 1, ndl
          ill = it + nstrtl - 1
          zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a2(it,il) * zbut(il)
        enddo
      enddo

    ! ltmtx == 2 means dense submatrix
    elseif (st_leafmtxp%st_lf(ip)%ltmtx == 2) then  ! (E1)
      do il = 1, ndl
        ill = il + nstrtl - 1
        do it = 1, ndt
          itt = it + nstrtt - 1
          zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a1(it,il) * zu(itt)
        enddo
      enddo
    endif

  ! time_matvec_loop_end = MPI_Wtime()
  ! time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start

  enddo

  ! print *, "Time of matvec loop:", time_matvec_loop

  do il = ls, le
  !$omp atomic
    zau(il) = zau(il) + zaut(il)
  enddo

end subroutine HACApK_adot_body_lfmtx_hyp_rise_64
! =============================================================================
! Rise 25 method stops here
! =============================================================================


! =============================================================================
! Rise's method 23 with inline comments on everything
! The goal of this method is to allow lower precision submatrices
! and check if the speed of matrix-vector computation increases.
! =============================================================================
! subroutine  HACApK_adot_body_lfmtx_hyp_rise_duplicate_cloning(zau, st_leafmtxp, st_ctl, zu, nd)
!   ! THIS FUNCTION IS VERY WRONG
!   ! The main subroutine that computes the matrix-vector multiplication
!   ! Input??:
!   !   zau: result of matvec -multiplication
!   !   st_leafmtxp: h-matrix
!   !   st_ctl: a lot of parameters
!   !   zu: source vector
!   !   nd: dimension of source vector
!   ! Returns??:
!   !   zau

!   implicit none
!   type(st_HACApK_leafmtxp) :: st_leafmtxp
!   type(st_HACApK_lcontrol) :: st_ctl
!   ! real(kind=8)
!   real*8 :: zau(nd), zu(nd)

!   ! Experiment 1 (E1): Changes zbut, zaut, a1 (and a2) into single precision
!   real*8, dimension(:), allocatable :: zbut
!   real*8, dimension(:), allocatable :: zaut

!   integer*4, pointer :: lpmd(:), lnp(:), lsp(:), ltmp(:)
!   1000 format(5(a, i10)/)
!   2000 format(5(a, f10.4)/)

!   integer, intent(in) :: nd

!   double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
!   integer :: mpinr, mpilog, nrank, icomm
!   integer :: ip, i, k, t,  il, ill, it, itt
!   integer :: nths, nthe, ith, ith1, ls, le
!   integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt
!   integer :: calcOrder

!   ! === Experiment 1 (Specification Section) ===
!   type(st_HACApK_leafmtxp_32) :: st_leafmtxp_32
!   real*4, dimension(:), allocatable :: zau_32, zu_32, zbut_32, zaut_32
!   !real*4, dimension(:, :), allocatable :: a1, a2
!   allocate(zau_32(nd))
!   allocate(zu_32(nd))
!   ! === end ===

!   lpmd => st_ctl%lpmd(:)
!   icomm = lpmd(1)
!   nrank = lpmd(2)
!   mpinr = lpmd(3)
!   mpilog = lpmd(4)
!   lnp(0:) => st_ctl%lnp
!   lsp(0:) => st_ctl%lsp
!   ltmp(0:) => st_ctl%lthr

!   nlf = st_leafmtxp%nlf
!   ktmax = st_leafmtxp%ktmax
!   ith = omp_get_thread_num()  ! thread_num
!   ith1 = ith+1
!   nths = ltmp(ith)
!   nthe = ltmp(ith1)-1
!   allocate(zaut(nd))
!   zaut(:) = 0.d0
!   allocate(zbut(ktmax))
!   ls = nd
!   le = 1

!   ! === Experiment 1 (Execution Section) ==
!   ! Experiment 1: Copies the whole 64b matrices into 32b versions and test speed
!   zau_32 = zau
!   zu_32 = zu
!   zbut_32 = zbut
!   zaut_32 = zaut


!   st_leafmtxp_32%nlf = st_leafmtxp%nlf
!   nlf = st_leafmtxp_32%nlf
!   st_leafmtxp_32%ktmax = st_leafmtxp%ktmax
!   ktmax = st_leafmtxp_32%ktmax
!   allocate(st_leafmtxp_32%st_lf(nlf))

!   do ip = nths, nthe

!     st_leafmtxp_32%st_lf(ip)%ltmtx = st_leafmtxp%st_lf(ip)%ltmtx
!     st_leafmtxp_32%st_lf(ip)%kt = st_leafmtxp%st_lf(ip)%kt
!     st_leafmtxp_32%st_lf(ip)%ndl = st_leafmtxp%st_lf(ip)%ndl
!     st_leafmtxp_32%st_lf(ip)%ndt = st_leafmtxp%st_lf(ip)%ndt
!     st_leafmtxp_32%st_lf(ip)%nstrtl = st_leafmtxp%st_lf(ip)%nstrtl
!     st_leafmtxp_32%st_lf(ip)%nstrtt = st_leafmtxp%st_lf(ip)%nstrtt


! !   do ip=1,nlf
!     ! -fbounds-check
!     ndl = st_leafmtxp_32%st_lf(ip)%ndl  ! (E1)
!     ndt = st_leafmtxp_32%st_lf(ip)%ndt  ! (E1)

! !  write(6,*) "ith=", ith, " ip=", ip

! !  stop

!     if (st_leafmtxp_32%st_lf(ip)%ltmtx == 1) then
!       kt = st_leafmtxp_32%st_lf(ip)%kt

! !      !$omp master

! !          write(6,*) "ip=", ip

!       allocate(st_leafmtxp_32%st_lf(ip)%a1(ndt, kt))
!       allocate(st_leafmtxp_32%st_lf(ip)%a2(ndl, kt))

! !      !$omp end master


!       do il = 1, kt
!         do it = 1 , ndt
!           st_leafmtxp_32%st_lf(ip)%a1(it,il) = st_leafmtxp%st_lf(ip)%a1(it,il)
!           ! st_leafmtxp_32%st_lf(ip)%a1(1:ndt, 1:kt) = st_leafmtxp%st_lf(ip)%a1(1:ndt, 1:kt)
!           ! test = st_leafmtxp%st_lf(ip)%a1(it, il)
!         enddo
!       enddo

!       do il = 1, kt
!         do it = 1 , ndl
!           st_leafmtxp_32%st_lf(ip)%a2(it,il) = st_leafmtxp%st_lf(ip)%a2(it,il)
!           ! st_leafmtxp_32%st_lf(ip)%a2(1:ndl, 1:kt) = st_leafmtxp%st_lf(ip)%a2(1:ndl, 1:kt)
!           ! test = st_leafmtxp%st_lf(ip)%a1(it, il)
!         enddo
!       enddo

!     elseif (st_leafmtxp_32%st_lf(ip)%ltmtx == 2) then
!       allocate(st_leafmtxp_32%st_lf(ip)%a1(ndt, ndl))
!       do il = 1, ndl
!         do it = 1 , ndt
!           st_leafmtxp_32%st_lf(ip)%a1(it,il) = st_leafmtxp%st_lf(ip)%a1(it,il)
!           ! st_leafmtxp_32%st_lf(ip)%a1(1:ndt,1:ndl) = st_leafmtxp%st_lf(ip)%a1(1:ndt,1:ndl)
!           ! test = st_leafmtxp%st_lf(ip)%a1(il, it)
!         enddo
!       enddo
!     endif
!   enddo

!   ! === semi ===

!   ! calcOrder = 0
!   time_matvec_loop = 0.0d0

!   do ip = nths, nthe
!     ! ndl = st_leafmtxp%st_lf(ip)%ndl  ! row size of submatrix
!     ! ndt = st_leafmtxp%st_lf(ip)%ndt  ! column size of submatrix
!     ! ns = ndl * ndt  ! This is never used??????
!     ! nstrtl = st_leafmtxp%st_lf(ip)%nstrtl  ! starting row no. of submatrix
!     ! nstrtt = st_leafmtxp%st_lf(ip)%nstrtt  ! starting col no. of submatrix

!     ndl = st_leafmtxp_32%st_lf(ip)%ndl  ! (E1)
!     ndt = st_leafmtxp_32%st_lf(ip)%ndt  ! (E1)
!     ns = ndl * ndt  ! (E1)
!     nstrtl = st_leafmtxp_32%st_lf(ip)%nstrtl  ! (E1)
!     nstrtt = st_leafmtxp_32%st_lf(ip)%nstrtt  ! (E1)

!     ! ! == E1 ==
!     ! if (st_leafmtxp_32%st_lf(ip)%ltmtx == 1) then  ! (E1)
!     !   allocate(a1_32(ndt, kt), a2_32(ndl, kt))
!     !   ! a1_32(1:ndt, 1:kt) = st_leafmtxp_32%st_lf(ip)%a1(1:ndt, 1:kt)
!     !   ! a2_32(1:ndl, 1:kt) = st_leafmtxp_32%st_lf(ip)%a2(1:ndl, 1:kt)
!     ! elseif (st_leafmtxp_32%st_lf(ip)%ltmtx == 2) then
!     !   allocate(a1_32(ndt, ndl))
!     !   ! a1_32(1:ndt, 1:ndl) = st_leafmtxp_32%st_lf(ip)%a1(1:ndt, 1:ndl)
!     ! endif
!     ! ! == == ==

!     ! We now have the coordinate of top left of submatrix (nstrtl & nstrtt)
!     ! We can now use that information and row & column size of submatrix to
!     ! loop over the whole thing and do something useful

!     if (nstrtl < ls) ls = nstrtl
!     if (nstrtl + ndl - 1 > le) le = nstrtl + ndl - 1

!     time_matvec_loop_start = MPI_Wtime()

!     ! ltmtx == 1 means low-rank submatrix
!     ! if (st_leafmtxp%st_lf(ip)%ltmtx == 1) then
!     if (st_leafmtxp_32%st_lf(ip)%ltmtx == 1) then  ! (E1)
!       ! Since we have decompose the submatrix into two matrices using
!       ! Adaptive Cross-Approximation (ACA) techniques, instead of just
!       ! size of (M X N), we now have two matrices of size (M X P) X (P X N)
!       ! Here specifically, (ndl X ndt) => (ndl X kt) X (kt X ndt)
!       ! kt = st_leafmtxp%st_lf(ip)%kt  ! (E1)
!       kt = st_leafmtxp_32%st_lf(ip)%kt  ! (E1)
!       ! zbut(1:kt) = 0.0d0  ! (Temporarily removed due to E1)
!       zbut_32(1:kt) = 0.0

!       ! We now do a double-decomposed-matrices-vector multiplication
!       ! (ndl X kt) X (kt X ndt) X (ndl)

!       ! The loop below is doing a (kt X ndt) X (ndl X 1) matvec_muliplication
!       ! The result is a (kt X 1) vector
!       do il = 1, kt
!         do it = 1 , ndt
!           itt = it + nstrtt - 1
!           ! zbut(il) = zbut(il) + st_leafmtxp%st_lf(ip)%a1(it,il) * zu(itt)  ! (rm by E1)
!           zbut_32(il) = zbut_32(il) + st_leafmtxp_32%st_lf(ip)%a1(it,il) * zu_32(itt)  ! (E1)
!           ! zbut_32(il) = zbut_32(il) + a1_32(it,il) * zu_32(itt)  ! (E1)
!         enddo
!       enddo

!       ! The loop below is doing a (ndl X kt) X (kt X 1) matvec_muliplication
!       ! The result is a (ndl X 1) vector
!       do il = 1, kt
!         do it = 1, ndl
!           ill = it + nstrtl - 1
!           ! zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a2(it,il) * zbut(il)  ! (rm by E1)
!           zaut_32(ill) = zaut_32(ill) + st_leafmtxp_32%st_lf(ip)%a2(it,il) * zbut_32(il)  ! (E1)
!           ! zaut_32(ill) = zaut_32(ill) + a2_32(it,il) * zbut_32(il)  ! (E1)
!         enddo
!       enddo

!     ! ltmtx == 2 means dense submatrix
!     ! elseif (st_leafmtxp%st_lf(ip)%ltmtx == 2) then
!     elseif (st_leafmtxp_32%st_lf(ip)%ltmtx == 2) then  ! (E1)
!       do il = 1, ndl
!         ill = il + nstrtl - 1
!         do it = 1, ndt
!           itt = it + nstrtt - 1
!           ! zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a1(it,il) * zu(itt)  ! (rm by E1)
!           zaut_32(ill) = zaut_32(ill) + st_leafmtxp_32%st_lf(ip)%a1(it,il) * zu_32(itt)  ! (E1)
!           ! zaut_32(ill) = zaut_32(ill) + a1_32(it,il) * zu_32(itt)  ! (E1)
!         enddo
!       enddo

!     endif

!     ! TODO: creates a ...%ltmtx == 5 single precision array

!     time_matvec_loop_end = MPI_Wtime()
!     time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start
!     ! calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
!   enddo

!   print *, "Time of matvec loop:", time_matvec_loop

!   do il = ls, le
!   !$omp atomic
!     zau(il) = zau(il) + zaut(il)  ! Not affected by E1
!   enddo

!   ! Question: What about deallocation??

! end subroutine HACApK_adot_body_lfmtx_hyp_rise_duplicate_cloning
! =============================================================================
! Rise 23 method stops here
! =============================================================================


!   method1
!***HACApK_adot_body_lfmtx_hyp_original
 subroutine HACApK_adot_body_lfmtx_hyp_original_ida(zau,st_leafmtxp,st_ctl,zu,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(*),zu(*)
 real*8,dimension(:),allocatable :: zbut
 real*8,dimension(:),allocatable :: zaut
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),ltmp(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

        integer, intent(in) :: nd

        double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
        integer :: mpinr, mpilog, nrank, icomm
        integer :: ip, i, k, t,  il, ill, it, itt
        integer :: nths, nthe, ith, ith1, ls, le
        integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt
!       integer :: calcOrder

!!!!!!!!!!!!!!!!!!!!!!!
!   Add by Iwashita
                double precision :: total_cost,workload,maxsize
                integer :: num_th
                integer :: kts1, kte1, kts2, kte2, thc, itemp, iii


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       calcOrder = 0
        time_matvec_loop = 0.0d0

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;ltmp(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 nlf=st_leafmtxp%nlf; ktmax=st_leafmtxp%ktmax
 ith = omp_get_thread_num()
 ith1 = ith+1
 nths=ltmp(ith); nthe=ltmp(ith1)-1
 allocate(zaut(nd)); zaut(:)=0.0d0
 allocate(zbut(ktmax))
 ls=nd; le=1

!!!!!!!!!!!!!!!!!!!!
! IWASHITA
 num_th=omp_get_num_threads()
! write(6,*)ltmp(0),ltmp(num_th)-1
 total_cost=0d0
do ip=ltmp(0),ltmp(num_th)-1
! do ip=nths,nthe
   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt

   if(nstrtl<ls) ls=nstrtl; if(nstrtl+ndl-1>le) le=nstrtl+ndl-1

   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt

!     zbut(1:kt)=0.0d0

!      total_cost=total_cost+kt

!     do il=1,kt
!       do it=1,ndt; itt=it+nstrtt-1
!         zbut(il)=zbut(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
!       enddo
!     enddo

      total_cost=total_cost+kt*ndt*2

!     do il=1,kt
!       do it=1,ndl; ill=it+nstrtl-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)
!       enddo
!     enddo

      total_cost=total_cost+kt*ndl*2

      if (maxsize<kt*(ndt+ndl)*2) then
        maxsize=kt*(ndt+ndl)*2
      endif


   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then

!     do il=1,ndl; ill=il+nstrtl-1
!       do it=1,ndt; itt=it+nstrtt-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
!       enddo
!     enddo

      total_cost=total_cost+ndt*ndl*2

   endif
!           time_matvec_loop_end = MPI_Wtime()
!           time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start
!           calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
 enddo

  if (ith==0) write(6,*) "total_cost=", total_cost
    workload=total_cost/num_th
  if (ith==0) then
     write(6,*) "workload=",workload
     write(6,*) "max size=",maxsize
  endif

thc=0
total_cost=0
itemp=0
kts1=-1
kte1=-1
kts2=-1
kte2=-1

do ip=ltmp(0),ltmp(num_th)-1

   if (thc==num_th-1) exit

   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt

   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     total_cost=total_cost+kt*ndt*2
     total_cost=total_cost+kt*ndl*2

   if (total_cost > workload) then

     itemp=((total_cost+1-workload)/(ndl+ndt)/2.0)

     if (itemp>0) then

       ltmp(thc+1)=ip
!       st_leafmtxp%st_lf(ip)%ltmtx=3

!     if (ith==thc) write(6,*) ith, total_cost, workload, ndl, ndt,itemp
!     write(6,*) "thc=",thc,"itemp=",itemp

       if (ith == thc) then
         kts2=1
         kte2=itemp
       endif

       if (ith == thc+1) then
         kts1=itemp+1
         kte1=kt
       endif

       total_cost=(kt-itemp)*(ndl+ndt)*2
       thc=thc+1

     else

       ltmp(thc+1)=ip+1
       total_cost=0
       thc=thc+1

     endif

   endif


   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then

      total_cost=total_cost+ndt*ndl*2

   if (total_cost > workload) then
     ltmp(thc+1)=ip+1
     total_cost=0
     thc=thc+1

   endif
   endif

enddo

!do ip=0, num_th-1
!$OMP BARRIER
!if (ith==ip) write(6,*) ith, ltmp(ith), ltmp(ith1)-1, kts1, kte1, kts2, kte2
!$OMP BARRIER
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!

if (ith.eq.0) time_matvec_start = MPI_Wtime()

do iii=1,100
 do ip=nths,nthe
   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
   if(nstrtl<ls) ls=nstrtl; if(nstrtl+ndl-1>le) le=nstrtl+ndl-1
!           time_matvec_loop_start = MPI_Wtime()
   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     zbut(1:kt)=0.0d0
     do il=1,kt
       do it=1,ndt; itt=it+nstrtt-1
         zbut(il)=zbut(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo
     do il=1,kt
       do it=1,ndl; ill=it+nstrtl-1
         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)
       enddo
     enddo
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo
   endif
!           time_matvec_loop_end = MPI_Wtime()
!           time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start
!           calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
 enddo
!       time_matvec_end = MPI_Wtime()
! deallocate(zbut)

!       time_reduction_start = MPI_Wtime()
 do il=ls,le
!$omp atomic
   zau(il)=zau(il)+zaut(il)
 enddo

!$OMP BARRIER
enddo

if (ith.eq.0)   time_reduction_end = MPI_Wtime()

if (ith.eq.0) write(6,*) "Elapsed",num_th,time_reduction_end-time_matvec_start

deallocate(zbut)


!       print '("myrank ", i3, "  time ", f)',  omp_get_thread_num(), time_matvec_end - time_matvec_start
!       print '("calcOrder: ", i0)', calcOrder


 end subroutine HACApK_adot_body_lfmtx_hyp_original_ida


!   method1
!***HACApK_adot_body_lfmtx_hyp_original
 subroutine HACApK_adot_body_lfmtx_hyp_original(zau,st_leafmtxp,st_ctl,zu,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(*),zu(*)
 real*8,dimension(:),allocatable :: zbut
 real*8,dimension(:),allocatable :: zaut
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),ltmp(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

        integer, intent(in) :: nd

        double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
        integer :: mpinr, mpilog, nrank, icomm
        integer :: ip, i, k, t,  il, ill, it, itt
        integer :: nths, nthe, ith, ith1, ls, le
        integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt
!       integer :: calcOrder

!!!!!!!!!!!!!!!!!!!!!!!
!   Add by Iwashita
                double precision :: total_cost,workload,maxsize
                integer :: num_th
                integer :: kts1, kte1, kts2, kte2, thc, itemp,iii


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       calcOrder = 0
        time_matvec_loop = 0.0d0

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;ltmp(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 nlf=st_leafmtxp%nlf; ktmax=st_leafmtxp%ktmax
 ith = omp_get_thread_num()
 ith1 = ith+1
 nths=ltmp(ith); nthe=ltmp(ith1)-1
 allocate(zaut(nd)); zaut(:)=0.0d0
 allocate(zbut(ktmax))
 ls=nd; le=1

!!!!!!!!!!!!!!!!!!!!
! IWASHITA
 num_th=omp_get_num_threads()
! write(6,*)ltmp(0),ltmp(num_th)-1
 total_cost=0d0
do ip=ltmp(0),ltmp(num_th)-1
! do ip=nths,nthe
   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt

   if(nstrtl<ls) ls=nstrtl; if(nstrtl+ndl-1>le) le=nstrtl+ndl-1

   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt

!     zbut(1:kt)=0.0d0

!      total_cost=total_cost+kt

!     do il=1,kt
!       do it=1,ndt; itt=it+nstrtt-1
!         zbut(il)=zbut(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
!       enddo
!     enddo

      total_cost=total_cost+kt*ndt*2

!     do il=1,kt
!       do it=1,ndl; ill=it+nstrtl-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)
!       enddo
!     enddo

      total_cost=total_cost+kt*ndl*2

      if (maxsize<kt*(ndt+ndl)*2) then
        maxsize=kt*(ndt+ndl)*2
      endif


   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then

!     do il=1,ndl; ill=il+nstrtl-1
!       do it=1,ndt; itt=it+nstrtt-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
!       enddo
!     enddo

      total_cost=total_cost+ndt*ndl*2

   endif
!           time_matvec_loop_end = MPI_Wtime()
!           time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start
!           calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
 enddo

!$OMP BARRIER

  if (ith==0) write(6,*) "total_cost=", total_cost
    workload=total_cost/num_th
  if (ith==0) then
     write(6,*) "workload=",workload
     write(6,*) "max size=",maxsize
  endif

thc=0
total_cost=0
itemp=0
kts1=-1
kte1=-1
kts2=-1
kte2=-1

do ip=ltmp(0),ltmp(num_th)-1
!$OMP BARRIER

   if (thc==num_th-1) exit

   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt

!   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
   if((st_leafmtxp%st_lf(ip)%ltmtx==1) .or. (st_leafmtxp%st_lf(ip)%ltmtx==3))then
     kt=st_leafmtxp%st_lf(ip)%kt
     total_cost=total_cost+kt*ndt*2
     total_cost=total_cost+kt*ndl*2


   if (total_cost > workload) then

!     write(6,*) ith,ip,"Kirikae"

     itemp=((total_cost+1-workload)/(ndl+ndt)/2.0)

     if (itemp>0) then

       ltmp(thc+1)=ip
       st_leafmtxp%st_lf(ip)%ltmtx=3

!     if (ith==thc) write(6,*) ith, total_cost, workload, ndl, ndt,itemp
!     write(6,*) "thc=",thc,"itemp=",itemp

!$OMP BARRIER

       if (ith == thc) then
         kts2=1
         kte2=itemp
!         write(6,*) "kts2,kte2", ith, kts2, kte2
       endif

       if (ith == thc+1) then
         kts1=itemp+1
         kte1=kt
!         write(6,*) "kts1,kte1", ith, kts1, kte1
       endif

       total_cost=(kt-itemp)*(ndl+ndt)*2
       thc=thc+1

     else

       ltmp(thc+1)=ip+1
       total_cost=0
       thc=thc+1

     endif

   endif


   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then

      total_cost=total_cost+ndt*ndl*2

   if (total_cost > workload) then
     ltmp(thc+1)=ip+1
     total_cost=0
     thc=thc+1

   endif
   endif

enddo
!$OMP BARRIER


!do ip=0, num_th-1
!!$OMP BARRIER
!if (ith==ip) write(6,*) ith, ltmp(ith), ltmp(ith1)-1, kts1, kte1, kts2, kte2
!!$OMP BARRIER
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!


ls=nd
le=1

 if (ith.eq.0) time_matvec_start = MPI_Wtime()
do iii=1,100


! do ip=nths,nthe
  do ip=ltmp(ith),ltmp(ith+1)-1
! IWASHITA
!!!!!

   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
   if(nstrtl<ls) ls=nstrtl; if(nstrtl+ndl-1>le) le=nstrtl+ndl-1
!           time_matvec_loop_start = MPI_Wtime()
   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     zbut(1:kt)=0.0d0
     do il=1,kt
       do it=1,ndt; itt=it+nstrtt-1
         zbut(il)=zbut(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo
     do il=1,kt
       do it=1,ndl; ill=it+nstrtl-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)
!$OMP ATOMIC
         zau(ill)=zau(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)

       enddo
     enddo
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
!$OMP ATOMIC
         zau(ill)=zau(ill)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)

       enddo
     enddo

!  IWASHITA
!
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==3)then

     kt=st_leafmtxp%st_lf(ip)%kt
     zbut(1:kt)=0.0d0

     do il=kts1,kte1
       do it=1,ndt; itt=it+nstrtt-1
         zbut(il)=zbut(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo

     do il=kts1,kte1
       do it=1,ndl; ill=it+nstrtl-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)
!$OMP ATOMIC
         zau(ill)=zau(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)


!         if (ith.eq.1) then
!           write(6,*) 'zaut',ill,zaut(ill),il
!         endif

       enddo
     enddo

!     write(6,*) "Enter 3",ith,ip,kts1,kte1,nstrtl,nstrtl+ndl-1


   endif
!           time_matvec_loop_end = MPI_Wtime()
!           time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start
!           calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
 enddo

   ip=ltmp(ith+1)
   if(st_leafmtxp%st_lf(ip)%ltmtx==3)then


     ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
     nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt

!     write(6,*) "thread number= ",ith,ip,"Last divided.",kts2,kte2,nstrtl,nstrtl+ndl-1

     if(nstrtl<ls) ls=nstrtl; if(nstrtl+ndl-1>le) le=nstrtl+ndl-1

     kt=st_leafmtxp%st_lf(ip)%kt
     zbut(1:kt)=0.0d0

     do il=kts2,kte2
       do it=1,ndt; itt=it+nstrtt-1
         zbut(il)=zbut(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo

     do il=kts2,kte2
       do it=1,ndl; ill=it+nstrtl-1
!         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)
!$OMP ATOMIC
zau(ill)=zau(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)
         if (ith.eq.0) then
!         write(6,*) 'zaut',ill,zaut(ill),il
         endif

       enddo
     enddo

   endif

!       time_matvec_end = MPI_Wtime()
! deallocate(zbut)

!       time_reduction_start = MPI_Wtime()
! do il=ls,le
!!$omp atomic
!   zau(il)=zau(il)+zaut(il)
! enddo

!$OMP Barrier
enddo

     if (ith.eq.0)  time_reduction_end = MPI_Wtime()

  if (ith.eq.0) write(6,*) "Elapsed time ",num_th,time_reduction_end-time_matvec_start

 deallocate(zbut)


!       print '("myrank ", i3, "  time ", f)',  omp_get_thread_num(), time_matvec_end - time_matvec_start
!       print '("calcOrder: ", i0)', calcOrder


 end subroutine HACApK_adot_body_lfmtx_hyp_original



!   method2
!***HACApK_adot_body_lfmtx_hyp
    subroutine HACApK_adot_body_lfmtx_hyp(zau, st_leafmtxp, st_ctl, zu, nd)
        type(st_HACApK_leafmtxp), intent(in) :: st_leafmtxp
        type(st_HACApK_lcontrol), intent(in) :: st_ctl
        real*8, intent(out) :: zau(*)
        real*8, intent(in)  :: zu(*)
        integer, intent(in) :: nd

        real*8, dimension(:), allocatable :: zaut, zbut
        real*8 :: time_start, time_bigmatrix_end, time_end, time_smallmatrix_start, time_smallmatrix_end
        integer*4, pointer :: ltmp(:)
        integer :: ip, k, t, it, itt, l, il, ill
!       integer :: calcOrder
        integer :: dividedBigmatrix
        integer :: ls, le, kt, ks, ke, nths, nthe
        integer :: ndl, ndt, nstrtl, nstrtt, ith, nThread, ipFirst, ipLast

        ith = omp_get_thread_num()
        nThread = omp_get_num_threads()
!       calcOrder = 0
        dividedBigmatrix = 0

        ltmp(0:) => st_ctl%lthr
        ipFirst = ltmp(0)
        ipLast  = ltmp(nThread) - 1

        allocate(zaut(nd), zbut(st_leafmtxp%ktmax))
        zaut(:) = 0.0d0
        ls = nd;    le = 1

        time_start = MPI_Wtime()

        do ip = ipFirst, ipLast
            ndl    = st_leafmtxp%st_lf(ip)%ndl;         ndt    = st_leafmtxp%st_lf(ip)%ndt
            nstrtl = st_leafmtxp%st_lf(ip)%nstrtl;      nstrtt = st_leafmtxp%st_lf(ip)%nstrtt
            kt  = st_leafmtxp%st_lf(ip)%kt

            if (st_leafmtxp%st_lf(ip)%ltmtx==3) then
                ks = loadBalanceQuot(kt, nThread) * ith + 1
                ke = loadBalanceQuot(kt, nThread) * (ith + 1)
                if(ke .gt. kt)      ke = kt
                zbut(1:kt) = 0.0d0

                if(nstrtl < ls)             ls = nstrtl
                if(nstrtl + ndl - 1 > le)   le = nstrtl + ndl - 1

                do k = ks, ke
                    do t = 1, ndt
                        itt = t + nstrtt - 1
                        zbut(k) = zbut(k) + st_leafmtxp%st_lf(ip)%a1(t, k) * zu(itt)
                    enddo
                enddo

                do k = ks, ke
                    do l = 1, ndl
                        ill = l + nstrtl - 1
                        zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a2(l, k) * zbut(k)
                    enddo
                enddo

!               dividedBigmatrix = dividedBigmatrix + 1
!               calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
            endif
        enddo

        time_bigmatrix_end = MPI_Wtime()
!!$omp barrier

        nths = ltmp(ith)
        nthe = ltmp(ith + 1) - 1

        time_smallmatrix_start = MPI_Wtime()
        do ip = nths, nthe
            ndl    = st_leafmtxp%st_lf(ip)%ndl;         ndt    = st_leafmtxp%st_lf(ip)%ndt
            nstrtl = st_leafmtxp%st_lf(ip)%nstrtl;      nstrtt = st_leafmtxp%st_lf(ip)%nstrtt
            kt  = st_leafmtxp%st_lf(ip)%kt

            if(nstrtl < ls)             ls = nstrtl
            if(nstrtl + ndl - 1 > le)   le = nstrtl + ndl - 1

            if(st_leafmtxp%st_lf(ip)%ltmtx==1) then
                zbut(1:kt) = 0.0d0
                do il = 1, kt
                    do it = 1, ndt
                        itt = it + nstrtt - 1
                        zbut(il) = zbut(il) + st_leafmtxp%st_lf(ip)%a1(it, il) * zu(itt)
                    enddo
                enddo
                do il = 1, kt
                    do it = 1, ndl
                        ill = it + nstrtl - 1
                        zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a2(it, il) * zbut(il)
                    enddo
                enddo
!               calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
            elseif(st_leafmtxp%st_lf(ip)%ltmtx==2) then
                do il = 1, ndl
                    ill = il + nstrtl - 1
                    do it = 1, ndt
                        itt = it + nstrtt - 1
                        zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a1(it, il) * zu(itt)
                    enddo
                enddo
!               calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
            endif
        enddo

        time_smallmatrix_end = MPI_Wtime()
        time_end = MPI_Wtime()


        do il = ls, le
!$omp atomic
            zau(il) = zau(il)+zaut(il)
        enddo

        print '("myrank ", i3, "  bigmatrix time ", f, "  smallmatrix time ", f, "  time ", f)',  omp_get_thread_num(),  time_bigmatrix_end - time_start, time_smallmatrix_end - time_smallmatrix_start, time_end - time_start
!       print '("myrank ", i3, "  divided bigmatrix: ", i0)',  omp_get_thread_num(),  dividedBigmatrix
!       print '("calcOrder: ", i0)', calcOrder

!$omp barrier
        deallocate(zaut, zbut)
    end subroutine HACApK_adot_body_lfmtx_hyp




!   method3
!***HACApK_adot_body_lfmtx_hyp_dynamic
    subroutine HACApK_adot_body_lfmtx_hyp_dynamic(zau, st_leafmtxp, st_ctl, zu, nd, chunksize)
        type(st_HACApK_leafmtxp), intent(in) :: st_leafmtxp
        type(st_HACApK_lcontrol), intent(in) :: st_ctl
        real*8, intent(out) :: zau(*)
        real*8, intent(in)  :: zu(*)
        integer, intent(in) :: nd, chunksize

        integer*4, pointer :: ltmp(:)
        integer :: ls, le, ip, ipFirst, ipLast, ndl, ndt, kt, nstrtl, nstrtt, il, ill, it, itt
!       integer :: calcOrder
        real*8, dimension(:), allocatable :: zaut, zbut

        ltmp(0:) => st_ctl%lthr
        ipFirst = ltmp(0)
        ipLast  = ltmp(omp_get_num_threads()) - 1

!       calcOrder = 0
        allocate(zaut(nd), zbut(st_leafmtxp%ktmax))
        zaut(:) = 0.0d0
        ls = nd;    le = 1

!$omp do schedule(dynamic, chunksize)
        do ip = ipFirst, ipLast
            ndl = st_leafmtxp%st_lf(ip)%ndl;            ndt = st_leafmtxp%st_lf(ip)%ndt
            nstrtl = st_leafmtxp%st_lf(ip)%nstrtl;      nstrtt = st_leafmtxp%st_lf(ip)%nstrtt
            kt  = st_leafmtxp%st_lf(ip)%kt

            if(nstrtl < ls)             ls = nstrtl
            if(nstrtl + ndl - 1 > le)   le = nstrtl + ndl - 1

            if(st_leafmtxp%st_lf(ip)%ltmtx==1) then
                zbut(1:kt) = 0.0d0
                do il = 1, kt
                    do it = 1, ndt
                        itt = it + nstrtt - 1
                        zbut(il) = zbut(il) + st_leafmtxp%st_lf(ip)%a1(it, il) * zu(itt)
                    enddo
                enddo
                do il = 1, kt
                    do it = 1, ndl
                        ill = it + nstrtl - 1
                        zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a2(it, il) * zbut(il)
                    enddo
                enddo
!               calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
            elseif(st_leafmtxp%st_lf(ip)%ltmtx==2) then
                do il = 1, ndl
                    ill = il + nstrtl - 1
                    do it = 1, ndt
                        itt = it + nstrtt - 1
                        zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a1(it, il) * zu(itt)
                    enddo
                enddo
!               calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
            endif
        enddo

        do il = ls, le
!$omp atomic
            zau(il) = zau(il)+zaut(il)
        enddo

!       print '("calcOrder: ", i0)', calcOrder

!$omp barrier
        deallocate(zaut, zbut)
    end subroutine HACApK_adot_body_lfmtx_hyp_dynamic



!   method4
!***HACApK_adot_body_lfmtx_hyp2
    subroutine HACApK_adot_body_lfmtx_hyp2(zau, st_leafmtxp, st_ctl, zu, nd, chunksize)
        type(st_HACApK_leafmtxp), intent(in) :: st_leafmtxp
        type(st_HACApK_lcontrol), intent(in) :: st_ctl
        real*8, intent(out) :: zau(*)
        real*8, intent(in)  :: zu(*)
        integer, intent(in) :: nd, chunksize

        real*8, dimension(:), allocatable :: zaut, zbut
        integer*4, pointer :: ltmp(:)
        integer :: ith, nThread, ip, ipFirst, ipLast
        integer :: ls, le, ks, ke, il, it, ill, itt, ndl, ndt, nstrtl, nstrtt, kt
!       integer :: calcOrder

        ith = omp_get_thread_num()
        nThread = omp_get_num_threads()
!       calcOrder = 0

        ltmp(0:) => st_ctl%lthr
        ipFirst = ltmp(0)
        ipLast  = ltmp(nThread) - 1

        allocate(zaut(nd), zbut(st_leafmtxp%ktmax))
        zaut(:) = 0.0d0
        ls = nd;    le = 1

        do ip = ipFirst, ipLast
            ndl    = st_leafmtxp%st_lf(ip)%ndl;         ndt    = st_leafmtxp%st_lf(ip)%ndt
            nstrtl = st_leafmtxp%st_lf(ip)%nstrtl;      nstrtt = st_leafmtxp%st_lf(ip)%nstrtt
            kt  = st_leafmtxp%st_lf(ip)%kt

            if(nstrtl < ls)             ls = nstrtl
            if(nstrtl + ndl - 1 > le)   le = nstrtl + ndl - 1

            if (st_leafmtxp%st_lf(ip)%ltmtx==3) then
                ks = loadBalanceQuot(kt, nThread) * ith + 1
                ke = loadBalanceQuot(kt, nThread) * (ith + 1)
                if(ke .gt. kt)      ke = kt
                zbut(1:kt) = 0.0d0

                do il = ks, ke
                    do it = 1, ndt
                        itt = it + nstrtt - 1
                        zbut(il) = zbut(il) + st_leafmtxp%st_lf(ip)%a1(it, il) * zu(itt)
                    enddo
                enddo

                do il = ks, ke
                    do it = 1, ndl
                        ill = it + nstrtl - 1
                        zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a2(it, il) * zbut(il)
                    enddo
                enddo

!               calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
            endif

        enddo
!$omp barrier

!$omp do schedule(dynamic, chunksize)
        do ip = ipFirst, ipLast
            ndl    = st_leafmtxp%st_lf(ip)%ndl;         ndt    = st_leafmtxp%st_lf(ip)%ndt
            nstrtl = st_leafmtxp%st_lf(ip)%nstrtl;      nstrtt = st_leafmtxp%st_lf(ip)%nstrtt
            kt  = st_leafmtxp%st_lf(ip)%kt

            if(nstrtl < ls)             ls = nstrtl
            if(nstrtl + ndl - 1 > le)   le = nstrtl + ndl - 1

            if(st_leafmtxp%st_lf(ip)%ltmtx==1) then
                zbut(1:kt) = 0.0d0
                do il = 1, kt
                    do it = 1, ndt
                        itt = it + nstrtt - 1
                        zbut(il) = zbut(il) + st_leafmtxp%st_lf(ip)%a1(it, il) * zu(itt)
                    enddo
                enddo
                do il = 1, kt
                    do it = 1, ndl
                        ill = it + nstrtl - 1
                        zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a2(it, il) * zbut(il)
                    enddo
                enddo

!               calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
            elseif(st_leafmtxp%st_lf(ip)%ltmtx==2) then
                do il = 1, ndl
                    ill = il + nstrtl - 1
                    do it = 1, ndt
                        itt = it + nstrtt - 1
                        zaut(ill) = zaut(ill) + st_leafmtxp%st_lf(ip)%a1(it, il) * zu(itt)
                    enddo
                enddo

!               calcOrder = calcOrder + leafcalcOrder(st_leafmtxp%st_lf(ip))
            endif
        enddo

        do il = ls, le
!$omp atomic
            zau(il) = zau(il)+zaut(il)
        enddo

!       print '("calcOrder: ", i0)', calcOrder

!$omp barrier
        deallocate(zaut, zbut)
    endsubroutine HACApK_adot_body_lfmtx_hyp2



!!!!!           st_lf(ndl, ndt)
!!!!!           a2 (ndl,kt) * a1 (ndt,kt)
!!!!!           ndt of a1 is col
!!!!!           ndl of a2 is row
!   method5
!***HACApK_adot_body_lfmtx_expand
    subroutine HACApK_adot_body_lfmtx_expand(zau, newLeafmtxp, zu, nd, chunksize)
        type(st_HACApK_leafmtxp), intent(in) :: newLeafmtxp
        real*8, intent(out) :: zau(*)
        real*8, intent(in)  :: zu(*)
        integer, intent(in) :: nd, chunksize

        !!      zbut:a1‚Æx‚ÌÏ      zaut:a2‚Æzbut‚ÌÏizau(=y)‚É‘«‚µ‚±‚Þj
        type(st_HACApK_leafmtx) :: newLeafmtx
        real*8, dimension(:), allocatable :: zaut, zbut
        double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end
        integer :: ip, i, k, t, l, il, ill, it, itt
        integer :: ndl, ndt, ns, kt, newKt, ktmax, nlf, newNlf, nstrtl, nstrtt
!       integer :: calcOrder
        integer :: ls, le, row, counter
        character :: filename * 30

        counter = 0
!       calcOrder = 0
        allocate(zaut(nd), zbut(newLeafmtxp%ktmax))
        zaut(1:nd) = 0.0d0
        zbut(1:newLeafmtxp%ktmax) = 0.0d0
        newNlf = newLeafmtxp%nlf
        ls = nd;    le = 1

        time_matvec_start = MPI_Wtime()
!$omp do schedule(dynamic, chunksize)
        do ip = 1, newNlf
            newLeafmtx = newLeafmtxp%st_lf(ip)
            ndl    = newLeafmtx%ndl;        ndt    = newLeafmtx%ndt
            nstrtl = newLeafmtx%nstrtl;     nstrtt = newLeafmtx%nstrtt
            newKt  = newLeafmtx%kt

            if(nstrtl < ls)             ls = nstrtl
            if(nstrtl + ndl - 1 > le)   le = nstrtl + ndl - 1

            if(newLeafmtx%ltmtx == 1) then
                zbut(1:newLeafmtxp%ktmax) = 0.0d0
                do k = 1, newKt
                    do t = 1, ndt
                        zbut(k) = zbut(k) + newLeafmtx%a1(t, k) * zu(nstrtt + t - 1)
                    enddo
                enddo
                do k = 1, newKt
                    do l = 1, ndl
                        zaut(nstrtl + l - 1) = zaut(nstrtl + l - 1) + newLeafmtx%a2(l, k) * zbut(k)
                    enddo
                enddo
            elseif(newLeafmtx%ltmtx == 2) then
                do l = 1, ndl
                    do t = 1, ndt
                        zaut(nstrtl + l - 1) = zaut(nstrtl + l - 1) + newLeafmtx%a1(t, l) * zu(nstrtt + t - 1)
                    enddo
                enddo
            endif

!           calcOrder = calcOrder + leafcalcOrder(newLeafmtx)
        enddo
        time_matvec_end = MPI_Wtime()

        time_reduction_start = MPI_Wtime()
        do il = ls, le
!$omp atomic
            zau(il) = zau(il) + zaut(il)
        enddo
        time_reduction_end = MPI_Wtime()

!       print '("calcOrder: ", i0)', calcOrder

        deallocate(zaut, zbut)
    end subroutine HACApK_adot_body_lfmtx_expand




!   method5B-2      methodNo:6
!***HACApK_adot_body_lfmtx_expandB2
    subroutine HACApK_adot_body_lfmtx_expandB2(zau, newLeafmtxp, zu, nd, chunksize, start, end)
        type(st_HACApK_leafmtxp), intent(in) :: newLeafmtxp
        real*8, intent(out) :: zau(*)
        real*8, intent(in)  :: zu(*)
        integer, intent(in) :: nd, chunksize, start, end

        type(st_HACApK_leafmtx) :: newLeafmtx
        real*8, dimension(:), allocatable :: zaut, zbut
        double precision :: time_matvec_start, time_matvec_end, time_reduction_start, time_reduction_end, time_matvec_loop_start, time_matvec_loop_end, time_matvec_loop
        integer :: ip, i, k, t, l, il, ill, it, itt
        integer :: ls, le, nths, nthe, ith, ith1
        integer :: ndl, ndt, ns, kt, ktmax, nlf, nstrtl, nstrtt, newKt
!       integer :: calcOrder

        allocate(zaut(nd), zbut(newLeafmtxp%ktmax))
        zaut(1:nd) = 0.0d0
        zbut(1:newLeafmtxp%ktmax) = 0.0d0
!       calcOrder = 0
        time_matvec_loop = 0.0d0
        ls = nd;    le = 1

!       time_matvec_start = MPI_Wtime()
        do ip = start, end
            newLeafmtx = newLeafmtxp%st_lf(ip)
            ndl    = newLeafmtx%ndl;        ndt    = newLeafmtx%ndt
            nstrtl = newLeafmtx%nstrtl;     nstrtt = newLeafmtx%nstrtt
            newKt  = newLeafmtx%kt

            if(nstrtl < ls)             ls = nstrtl
            if(nstrtl + ndl - 1 > le)   le = nstrtl + ndl - 1

!           time_matvec_loop_start = MPI_Wtime()
            if(newLeafmtx%ltmtx == 1) then
                zbut(1:newLeafmtxp%ktmax) = 0.0d0
                do k = 1, newKt
                    do t = 1, ndt
                        zbut(k) = zbut(k) + newLeafmtx%a1(t, k) * zu(nstrtt + t - 1)
                    enddo
                enddo
                do k = 1, newKt
                    do l = 1, ndl
                        zaut(nstrtl + l - 1) = zaut(nstrtl + l - 1) + newLeafmtx%a2(l, k) * zbut(k)
                    enddo
                enddo
            elseif(newLeafmtx%ltmtx == 2) then
                do l = 1, ndl
                    do t = 1, ndt
                        zaut(nstrtl + l - 1) = zaut(nstrtl + l - 1) + newLeafmtx%a1(t, l) * zu(nstrtt + t - 1)
                    enddo
                enddo
            endif
!           time_matvec_loop_end = MPI_Wtime()
!           time_matvec_loop = time_matvec_loop + time_matvec_loop_end - time_matvec_loop_start

!           calcOrder = calcOrder + leafcalcOrder(newLeafmtx)
        enddo
!       time_matvec_end = MPI_Wtime()
!$omp barrier

!       time_reduction_start = MPI_Wtime()
        do il = ls, le
!$omp atomic
            zau(il) = zau(il) + zaut(il)
        enddo
!       time_reduction_end = MPI_Wtime()

!       print '("calcOrder: ", i0)', calcOrder



!$omp barrier
        deallocate(zaut, zbut)
    end subroutine HACApK_adot_body_lfmtx_expandB2





end module m_HACApK_solve_methods
