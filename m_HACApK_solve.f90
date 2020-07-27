!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 1.0.0                                             *
!                                                                     *
!   License                                                           *
!     This file is part of HACApK.                                    *
!     HACApK is a free software, you can use it under the terms       *
!     of The MIT License (MIT). See LICENSE file and User's guide     *
!     for more details.                                               *
!                                                                     *
!   ppOpen-HPC project:                                               *
!     Open Source Infrastructure for Development and Execution of     *
!     Large-Scale Scientific Applications on Post-Peta-Scale          *
!     Supercomputers with Automatic Tuning (AT).                      *
!                                                                     *
!   Sponsorship:                                                      *
!     Japan Science and Technology Agency (JST), Basic Research       *
!     Programs: CREST, Development of System Software Technologies    *
!     for post-Peta Scale High Performance Computing.                 *
!                                                                     *
!   Copyright (c) 2015 <Akihiro Ida and Takeshi Iwashita>             *
!                                                                     *
!=====================================================================*
!C***********************************************************************
!C  This file includes routines for utilizing H-matrices, such as solving
!C  linear system with an H-matrix as the coefficient matrix and
!C  multiplying an H-matrix and a vector,
!C  created by Akihiro Ida at Kyoto University on May 2012
!C***********************************************************************
!C***********************************************************************
!C  Creating a new branch under Rise in 2018
!C***********************************************************************
module m_HACApK_solve
        use m_HACApK_base
        use m_HACApK_mine
        use m_HACApK_solve_methods
        implicit real*8(a-h, o-z)
        implicit integer*4(i-n)

contains


!***HACApK_bicgstab_lfmtx
    subroutine HACApK_bicgstab_lfmtx(st_leafmtxp, st_ctl, u, b, param, nd, nstp, lrtrn)
        include 'mpif.h'
        type(st_HACApK_leafmtxp), intent(in) :: st_leafmtxp
        type(st_HACApK_lcontrol), intent(in) :: st_ctl
        real*8 :: u(nd), b(nd)
        real*8, intent(in) :: param(*)

        real*8, dimension(:), allocatable :: zr, zshdw, zp, zt, zkp, zakp, zkt, zakt
        integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
        1000 format(5(a, i10)/)
        2000 format(5(a, f10.4)/)

!       integer*4 :: mpinr, icomm, ierr, mstep, in
!       real*8 :: st_measure_time, en_measure_time, time
        real*8 :: eps, bnorm, zrnorm, znom, znomold, alpha=0.0, beta=0.0, zeta=0.0

        lpmd => st_ctl%lpmd(:)
        mpinr = lpmd(3)
        icomm = lpmd(1)

        call MPI_Barrier(icomm, ierr)
        st_measure_time = MPI_Wtime()
        if(st_ctl%param(1)>0 .and. mpinr==0)        print *, 'HACApK_bicgstab_lfmtx start'

        mstep = param(83)
        eps = param(91)

        allocate(zr(nd), zshdw(nd), zp(nd), zt(nd), zkp(nd), zakp(nd), zkt(nd), zakt(nd))
        zp(1:nd) = 0.0d0
        zakp(1:nd) = 0.0d0
        bnorm = dsqrt(HACApK_dotp_d(nd, b, b))
        zr(:nd) = b(:nd)

        call HACApK_adotsub_lfmtx_p(zr, st_leafmtxp, st_ctl, u, nd)
        zshdw(:nd) = zr(:nd)
        zrnorm = dsqrt(HACApK_dotp_d(nd, zr, zr))
        if(st_ctl%param(1)>0 .and. mpinr==0)        print *, 'Original relative residual norm =', zrnorm/bnorm
        if(zrnorm/bnorm < eps)      return

! mstep = 1
        do in = 1, mstep
            zp(:nd) = zr(:nd) + beta * (zp(:nd) - zeta * zakp(:nd))
            zkp(:nd) = zp(:nd)
            call HACApK_adot_lfmtx_p(zakp, st_leafmtxp, st_ctl, zkp, nd)
! exit
            znom = HACApK_dotp_d(nd, zshdw, zr)
            alpha = znom / HACApK_dotp_d(nd, zshdw, zakp)           ! znom/zden
            znomold = znom
            zt(:nd) = zr(:nd) - alpha * zakp(:nd)
            zkt(:nd) = zt(:nd)
            call HACApK_adot_lfmtx_p(zakt, st_leafmtxp, st_ctl, zkt, nd)

            znom = HACApK_dotp_d(nd, zakt, zt)
            zeta = znom / HACApK_dotp_d(nd, zakt, zakt)             ! znom/zden
            u(:nd) = u(:nd) + alpha * zkp(:nd) + zeta * zkt(:nd)
            zr(:nd) = zt(:nd) - zeta * zakt(:nd)
            beta = alpha/zeta * HACApK_dotp_d(nd, zshdw, zr) / znomold
            zrnorm = dsqrt(HACApK_dotp_d(nd, zr, zr))
            call MPI_Barrier(icomm, ierr)

            en_measure_time = MPI_Wtime()
            time = en_measure_time - st_measure_time
            if(st_ctl%param(1)>0 .and. mpinr==0)        print *, in, time, log10(zrnorm/bnorm)
            if(zrnorm/bnorm < eps)      exit
        enddo
    end subroutine HACApK_bicgstab_lfmtx


!***HACApK_adotsub_lfmtx_p
    subroutine HACApK_adotsub_lfmtx_p(zr, st_leafmtxp, st_ctl, zu, nd)
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zu(nd), zr(nd)
        real*8, dimension(:), allocatable :: zau
        integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)

        allocate(zau(nd))

        call HACApK_adot_lfmtx_p(zau, st_leafmtxp, st_ctl, zu, nd)
        zr(1:nd) = zr(1:nd) - zau(1:nd)
        deallocate(zau)
    end subroutine HACApK_adotsub_lfmtx_p


!***HACApK_adot_lfmtx_p
    subroutine HACApK_adot_lfmtx_p(zau, st_leafmtxp, st_ctl, zu, nd)
        include 'mpif.h'
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zau(nd), zu(nd)
        real*8, dimension(:), allocatable :: wws, wwr
        integer*4 :: ISTATUS(MPI_STATUS_SIZE), isct(2), irct(2)
        integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
        1000 format(5(a, i10)/)
        2000 format(5(a, f10.4)/)

        lpmd => st_ctl%lpmd(:)
        lnp(0:) => st_ctl%lnp
        lsp(0:) => st_ctl%lsp
        mpinr = lpmd(3)
        mpilog = lpmd(4)
        nrank = lpmd(2)
        icomm = lpmd(1)
        allocate(wws(maxval(lnp(0:nrank-1))), wwr(maxval(lnp(0:nrank-1))))
        zau(:) = 0.0d0

        call HACApK_adot_body_lfmtx(zau, st_leafmtxp, st_ctl, zu, nd)
        if(nrank==1) return

        wws(1:lnp(mpinr)) = zau(lsp(mpinr):lsp(mpinr) + lnp(mpinr)-1)
        ncdp = mod(mpinr + 1, nrank)
        ncsp = mod(mpinr + nrank-1, nrank)
! write(mpilog, 1000) 'destination process = ', ncdp, '; source process = ', ncsp
        isct(1) = lnp(mpinr)
        isct(2) = lsp(mpinr)
! irct = lnp(ncsp)
        do ic = 1, nrank-1
!   idp = mod(mpinr + ic, nrank) ! rank of destination process
!   isp = mod(mpinr + nrank + ic-2, nrank) ! rank of source process
            call MPI_SENDRECV(isct, 2, MPI_INTEGER, ncdp, 1, irct, 2, MPI_INTEGER, ncsp, 1, icomm, ISTATUS, ierr)
!   write(mpilog, 1000) 'ISTATUS = ', ISTATUS, '; ierr = ', ierr
!   write(mpilog, 1000) 'ic = ', ic, '; isct = ', isct(1), '; irct = ', irct(1), '; ivsps = ', isct(2), '; ivspr = ', irct(2)
            call MPI_SENDRECV(wws, isct, MPI_DOUBLE_PRECISION, ncdp, 1, wwr, irct, MPI_DOUBLE_PRECISION, ncsp, 1,icomm, ISTATUS, ierr)

!   write(mpilog, 1000) 'ISTATUS = ', ISTATUS, '; ierr = ', ierr
            zau(irct(2):irct(2) + irct(1)-1) = zau(irct(2):irct(2) + irct(1)-1) + wwr(:irct(1))
            wws(:irct(1)) = wwr(:irct(1))
            isct = irct
!   write(mpilog, 1000) 'ic = ', ic, '; isct = ', isct
        enddo

        deallocate(wws, wwr)
    end subroutine HACApK_adot_lfmtx_p


!***HACApK_adot_body_lfmtx
    RECURSIVE subroutine HACApK_adot_body_lfmtx(zau, st_leafmtxp, st_ctl, zu, nd)
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zau(nd), zu(nd)
        real*8, dimension(:), allocatable :: zbu

        nlf = st_leafmtxp%nlf

        do ip = 1, nlf
            ndl = st_leafmtxp%st_lf(ip)%ndl
            ndt = st_leafmtxp%st_lf(ip)%ndt
            nstrtl = st_leafmtxp%st_lf(ip)%nstrtl
            nstrtt = st_leafmtxp%st_lf(ip)%nstrtt

            if(st_leafmtxp%st_lf(ip)%ltmtx==1) then
                kt = st_leafmtxp%st_lf(ip)%kt
                allocate(zbu(kt))
                zbu(:) = 0.0d0
                do il = 1, kt
                    do it = 1, ndt
                        itt = it + nstrtt-1
                        zbu(il) = zbu(il) + st_leafmtxp%st_lf(ip)%a1(it, il) * zu(itt)
                    enddo
                enddo
                do il = 1, kt
                    do it = 1, ndl
                        ill = it + nstrtl-1
                        zau(ill) = zau(ill) + st_leafmtxp%st_lf(ip)%a2(it, il) * zbu(il)
                    enddo
                enddo
                deallocate(zbu)
            elseif(st_leafmtxp%st_lf(ip)%ltmtx==2) then
                do il = 1, ndl
                    ill = il + nstrtl-1
                    do it = 1, ndt
                        itt = it + nstrtt-1
                        zau(ill) = zau(ill) + st_leafmtxp%st_lf(ip)%a1(it, il) * zu(itt)
                    enddo
                enddo
            endif
        enddo
    end subroutine HACApK_adot_body_lfmtx

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!***HACApK_bicgstab_lfmtx_hyp
    subroutine HACApK_bicgstab_lfmtx_hyp(st_leafmtxp, st_ctl, u, b, param, nd, nstp, lrtrn)
        include 'mpif.h'
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: u(nd), b(nd)
        real*8 :: param(*)

        real*8, dimension(:), allocatable :: zr, zshdw, zp, zt, zkp, zakp, zkt, zakt
        real*8, dimension(:), allocatable :: wws, wwr
        integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
        integer*4 :: isct(2), irct(2)
        integer :: mstep, in
        real*8  :: alpha=0.0, beta=0.0, zeta=0.0
        real*8, save :: time_total=0.0


        ! RISE PART STARTS HERE
        ! type(st_HACApK_leafmtxp_32) :: st_leafmtxp_32
        type(st_HACApK_leafmtxp_scaled) :: st_leafmtxp_scaled
        ! type(st_HACApK_leafmtxp_scaled_32) :: st_leafmtxp_scaled_32
        type(st_HACApK_leafmtxp_mixed) :: st_leafmtxp_mixed

        real*4 :: placeholder_32(nd)
        real*8, dimension(nd) :: dummy_64
        real*4, dimension(nd) :: dummy_32
        real*8 :: max, SPLIT_CRITERION
        integer*4 :: split_value, track_64, track_32

        dummy_64(:nd) = 1d-6
        dummy_32(:nd) = 1e-6

        ! =====================================================================
        ! =====================================================================
        ! =====================================================================
        ! TODO:
        ! 1. Rise: Low-rank -> 32b, Full-rank -> 64b (Phase 3 - b??)
        ! 2. Iwashita: Low-rank, according to diag, precision changes. -> simplify (removal), 32b, and 64b
        ! 3. Iwashita: Full-rank, how to change precision?? Dunno, ask sensei
        ! Another idea: Full-rank -> 32b. Low-rank 32-16?
        ! Can we even decide which lesser full-ranks to be 32b and which full-ranks to be 64b?

        ! Thesis TODO:
        ! 1. Understanding Phase 2: Turn off O3 (O0) in Makefile and see the time-threadcount behavior
        ! 2. Speeding Up Phase 2: Allocate Diagonal in the middle of a1 and a2, not after
        ! =====================================================================
        ! =====================================================================
        ! =====================================================================

        ! ! copy 64 bit into 32 bit
        ! st_leafmtxp_32%nd = st_leafmtxp%nd
        ! st_leafmtxp_32%nlf = st_leafmtxp%nlf
        ! st_leafmtxp_32%nlfkt = st_leafmtxp%nlfkt
        ! st_leafmtxp_32%ktmax = st_leafmtxp%ktmax
        ! allocate(st_leafmtxp_32%st_lf(st_leafmtxp_32%nlf))

        ! do ip = 1, st_leafmtxp_32%nlf
        !     st_leafmtxp_32%st_lf(ip)%ltmtx = st_leafmtxp%st_lf(ip)%ltmtx
        !     st_leafmtxp_32%st_lf(ip)%ndl = st_leafmtxp%st_lf(ip)%ndl
        !     st_leafmtxp_32%st_lf(ip)%ndt = st_leafmtxp%st_lf(ip)%ndt
        !     st_leafmtxp_32%st_lf(ip)%nstrtl = st_leafmtxp%st_lf(ip)%nstrtl
        !     st_leafmtxp_32%st_lf(ip)%nstrtt = st_leafmtxp%st_lf(ip)%nstrtt

        !     if (st_leafmtxp_32%st_lf(ip)%ltmtx == 1) then
        !         st_leafmtxp_32%st_lf(ip)%kt = st_leafmtxp%st_lf(ip)%kt

        !         allocate(st_leafmtxp_32%st_lf(ip)%a1(st_leafmtxp_32%st_lf(ip)%ndt, st_leafmtxp_32%st_lf(ip)%kt))
        !         allocate(st_leafmtxp_32%st_lf(ip)%a2(st_leafmtxp_32%st_lf(ip)%ndl, st_leafmtxp_32%st_lf(ip)%kt))
        !         st_leafmtxp_32%st_lf(ip)%a1(1:st_leafmtxp_32%st_lf(ip)%ndt, 1:st_leafmtxp_32%st_lf(ip)%kt) = st_leafmtxp%st_lf(ip)%a1(1:st_leafmtxp%st_lf(ip)%ndt, 1:st_leafmtxp%st_lf(ip)%kt)
        !         st_leafmtxp_32%st_lf(ip)%a2(1:st_leafmtxp_32%st_lf(ip)%ndl, 1:st_leafmtxp_32%st_lf(ip)%kt) = st_leafmtxp%st_lf(ip)%a2(1:st_leafmtxp%st_lf(ip)%ndl, 1:st_leafmtxp%st_lf(ip)%kt)

        !     elseif (st_leafmtxp_32%st_lf(ip)%ltmtx == 2) then
        !         allocate(st_leafmtxp_32%st_lf(ip)%a1(st_leafmtxp_32%st_lf(ip)%ndt, st_leafmtxp_32%st_lf(ip)%ndl))
        !         st_leafmtxp_32%st_lf(ip)%a1(1:st_leafmtxp_32%st_lf(ip)%ndt, 1:st_leafmtxp_32%st_lf(ip)%ndl) = st_leafmtxp%st_lf(ip)%a1(1:st_leafmtxp%st_lf(ip)%ndt,1:st_leafmtxp%st_lf(ip)%ndl)

        !     endif

        ! enddo

        ! ! copy original to the modified one that includes an additional diagonal array (for scaling)
        st_leafmtxp_scaled%nd = st_leafmtxp%nd
        st_leafmtxp_scaled%nlf = st_leafmtxp%nlf
        st_leafmtxp_scaled%nlfkt = st_leafmtxp%nlfkt
        st_leafmtxp_scaled%ktmax = st_leafmtxp%ktmax
        allocate(st_leafmtxp_scaled%st_lf(st_leafmtxp_scaled%nlf))

        do ip = 1, st_leafmtxp_scaled%nlf
            st_leafmtxp_scaled%st_lf(ip)%ltmtx = st_leafmtxp%st_lf(ip)%ltmtx
            st_leafmtxp_scaled%st_lf(ip)%ndl = st_leafmtxp%st_lf(ip)%ndl
            st_leafmtxp_scaled%st_lf(ip)%ndt = st_leafmtxp%st_lf(ip)%ndt
            st_leafmtxp_scaled%st_lf(ip)%nstrtl = st_leafmtxp%st_lf(ip)%nstrtl
            st_leafmtxp_scaled%st_lf(ip)%nstrtt = st_leafmtxp%st_lf(ip)%nstrtt

            if (st_leafmtxp_scaled%st_lf(ip)%ltmtx == 1) then
                st_leafmtxp_scaled%st_lf(ip)%kt = st_leafmtxp%st_lf(ip)%kt

                allocate(st_leafmtxp_scaled%st_lf(ip)%a1(st_leafmtxp_scaled%st_lf(ip)%ndt, st_leafmtxp_scaled%st_lf(ip)%kt))
                allocate(st_leafmtxp_scaled%st_lf(ip)%a2(st_leafmtxp_scaled%st_lf(ip)%ndl, st_leafmtxp_scaled%st_lf(ip)%kt))
                st_leafmtxp_scaled%st_lf(ip)%a1(1:st_leafmtxp_scaled%st_lf(ip)%ndt, 1:st_leafmtxp_scaled%st_lf(ip)%kt) = st_leafmtxp%st_lf(ip)%a1(1:st_leafmtxp%st_lf(ip)%ndt, 1:st_leafmtxp%st_lf(ip)%kt)
                st_leafmtxp_scaled%st_lf(ip)%a2(1:st_leafmtxp_scaled%st_lf(ip)%ndl, 1:st_leafmtxp_scaled%st_lf(ip)%kt) = st_leafmtxp%st_lf(ip)%a2(1:st_leafmtxp%st_lf(ip)%ndl, 1:st_leafmtxp%st_lf(ip)%kt)

            elseif (st_leafmtxp_scaled%st_lf(ip)%ltmtx == 2) then
                allocate(st_leafmtxp_scaled%st_lf(ip)%a1(st_leafmtxp_scaled%st_lf(ip)%ndt, st_leafmtxp_scaled%st_lf(ip)%ndl))
                st_leafmtxp_scaled%st_lf(ip)%a1(1:st_leafmtxp_scaled%st_lf(ip)%ndt, 1:st_leafmtxp_scaled%st_lf(ip)%ndl) = st_leafmtxp%st_lf(ip)%a1(1:st_leafmtxp%st_lf(ip)%ndt,1:st_leafmtxp%st_lf(ip)%ndl)

            endif

        enddo


        ! allocate diagonal and scale the st_leafmtxp_scaled
        do ip = 1, st_leafmtxp_scaled%nlf
            if (st_leafmtxp_scaled%st_lf(ip)%ltmtx == 1) then
                allocate(st_leafmtxp_scaled%st_lf(ip)%diag(st_leafmtxp_scaled%st_lf(ip)%kt))
                st_leafmtxp_scaled%st_lf(ip)%diag(1:st_leafmtxp_scaled%st_lf(ip)%kt) = 1d0

                ! It is possible to merge the following two nested loops together
                ! It is also possible to parallelize the operations
                ! But since we are not timing these, no optimizations needed for now
                ! In the future, we will do these straight away in the generation phase

                ! find max of a1's columns, scale it, and add to diag
                do k = 1, st_leafmtxp_scaled%st_lf(ip)%kt
                    max = 0

                    ! find max of each column part
                    do n = 1, st_leafmtxp_scaled%st_lf(ip)%ndt
                        if (abs(st_leafmtxp_scaled%st_lf(ip)%a1(n, k)) > max) then
                            max = abs(st_leafmtxp_scaled%st_lf(ip)%a1(n, k))
                        endif
                    enddo

                    ! actually the following if/else is not needed because useless zero rows/columns do not exist
                    if (max /= 0) then
                        ! scale each column part
                        do n = 1, st_leafmtxp_scaled%st_lf(ip)%ndt
                            st_leafmtxp_scaled%st_lf(ip)%a1(n, k) = 1.0 * st_leafmtxp_scaled%st_lf(ip)%a1(n, k) / max
                        enddo
                        ! add to diag part
                        st_leafmtxp_scaled%st_lf(ip)%diag(k) = st_leafmtxp_scaled%st_lf(ip)%diag(k) * max
                    ! else
                        ! st_leafmtxp_scaled%st_lf(ip)%diag(k) = 0d0
                    endif
                enddo

                ! repeat for a2: find max of a2's rows, scale it, and add to diag
                do k = 1, st_leafmtxp_scaled%st_lf(ip)%kt
                    max = 0
                    ! find max of each row part
                    do n = 1, st_leafmtxp_scaled%st_lf(ip)%ndl
                        if (abs(st_leafmtxp_scaled%st_lf(ip)%a2(n, k)) > max) then
                            max = abs(st_leafmtxp_scaled%st_lf(ip)%a2(n, k))
                        endif
                    enddo

                    ! actually the following if/else is not needed because useless zero rows/columns do not exist
                    if (max /= 0) then
                        ! scale each row part
                        do n = 1, st_leafmtxp_scaled%st_lf(ip)%ndl
                            st_leafmtxp_scaled%st_lf(ip)%a2(n, k) = 1.0 * st_leafmtxp_scaled%st_lf(ip)%a2(n, k) / max
                        enddo
                        ! add to diag part
                        st_leafmtxp_scaled%st_lf(ip)%diag(k) = st_leafmtxp_scaled%st_lf(ip)%diag(k) * max
                    ! else
                        ! st_leafmtxp_scaled%st_lf(ip)%diag(k) = 0d0
                    endif
                enddo

            endif
        enddo


        ! ! copy scaled 64 bit into scaled 32 bit
        ! st_leafmtxp_scaled_32%nd = st_leafmtxp_scaled%nd
        ! st_leafmtxp_scaled_32%nlf = st_leafmtxp_scaled%nlf
        ! st_leafmtxp_scaled_32%nlfkt = st_leafmtxp_scaled%nlfkt
        ! st_leafmtxp_scaled_32%ktmax = st_leafmtxp_scaled%ktmax
        ! allocate(st_leafmtxp_scaled_32%st_lf(st_leafmtxp_scaled_32%nlf))

        ! do ip = 1, st_leafmtxp_scaled_32%nlf
        !     st_leafmtxp_scaled_32%st_lf(ip)%ltmtx = st_leafmtxp_scaled%st_lf(ip)%ltmtx
        !     st_leafmtxp_scaled_32%st_lf(ip)%ndl = st_leafmtxp_scaled%st_lf(ip)%ndl
        !     st_leafmtxp_scaled_32%st_lf(ip)%ndt = st_leafmtxp_scaled%st_lf(ip)%ndt
        !     st_leafmtxp_scaled_32%st_lf(ip)%nstrtl = st_leafmtxp_scaled%st_lf(ip)%nstrtl
        !     st_leafmtxp_scaled_32%st_lf(ip)%nstrtt = st_leafmtxp_scaled%st_lf(ip)%nstrtt

        !     if (st_leafmtxp_scaled_32%st_lf(ip)%ltmtx == 1) then
        !         st_leafmtxp_scaled_32%st_lf(ip)%kt = st_leafmtxp_scaled%st_lf(ip)%kt

        !         allocate(st_leafmtxp_scaled_32%st_lf(ip)%a1(st_leafmtxp_scaled_32%st_lf(ip)%ndt, st_leafmtxp_scaled_32%st_lf(ip)%kt))
        !         allocate(st_leafmtxp_scaled_32%st_lf(ip)%a2(st_leafmtxp_scaled_32%st_lf(ip)%ndl, st_leafmtxp_scaled_32%st_lf(ip)%kt))
        !         allocate(st_leafmtxp_scaled_32%st_lf(ip)%diag(st_leafmtxp_scaled_32%st_lf(ip)%kt))

        !         st_leafmtxp_scaled_32%st_lf(ip)%a1(1:st_leafmtxp_scaled_32%st_lf(ip)%ndt, 1:st_leafmtxp_scaled_32%st_lf(ip)%kt) = st_leafmtxp_scaled%st_lf(ip)%a1(1:st_leafmtxp_scaled%st_lf(ip)%ndt, 1:st_leafmtxp_scaled%st_lf(ip)%kt)
        !         st_leafmtxp_scaled_32%st_lf(ip)%a2(1:st_leafmtxp_scaled_32%st_lf(ip)%ndl, 1:st_leafmtxp_scaled_32%st_lf(ip)%kt) = st_leafmtxp_scaled%st_lf(ip)%a2(1:st_leafmtxp_scaled%st_lf(ip)%ndl, 1:st_leafmtxp_scaled%st_lf(ip)%kt)
        !         st_leafmtxp_scaled_32%st_lf(ip)%diag(1:st_leafmtxp_scaled_32%st_lf(ip)%kt) = st_leafmtxp_scaled%st_lf(ip)%diag(1:st_leafmtxp_scaled%st_lf(ip)%kt)

        !     elseif (st_leafmtxp_scaled_32%st_lf(ip)%ltmtx == 2) then

        !         allocate(st_leafmtxp_scaled_32%st_lf(ip)%a1(st_leafmtxp_scaled_32%st_lf(ip)%ndt, st_leafmtxp_scaled_32%st_lf(ip)%ndl))
        !         st_leafmtxp_scaled_32%st_lf(ip)%a1(1:st_leafmtxp_scaled_32%st_lf(ip)%ndt, 1:st_leafmtxp_scaled_32%st_lf(ip)%ndl) = st_leafmtxp_scaled%st_lf(ip)%a1(1:st_leafmtxp_scaled%st_lf(ip)%ndt,1:st_leafmtxp_scaled%st_lf(ip)%ndl)

        !     endif

        ! enddo

        ! ============================================== START OF PHASE 3 ================================================
        ! NEED TO CHECK THE SMALLEST DIAG, IF LESS THAN 1E-6 EXIST THEN MAYBE WE SHOULD DO SIMPLIFICATION, IF NOT, DON'T WASTE TIME

        ! 0. Analysis Check
        ! !$omp master
        ! ! print *, 'DEBUG 64'

        ! open(10, file='analysis_check.txt')
        !   do ip = 1, st_leafmtxp_scaled%nlf
        !     if (st_leafmtxp_scaled%st_lf(ip)%ltmtx == 1) then
        !         do k = 1, st_leafmtxp_scaled%st_lf(ip)%kt
        !             write(10, *) st_leafmtxp_scaled%st_lf(ip)%diag(k)
        !         enddo
        !     endif
        !   enddo
        ! close(10)

        ! !$omp end master
        ! !$omp barrier

        ! stop

        ! 1. Scale original, get diagonal
        ! For our convenience, we use scaled instead of redoing
        st_leafmtxp_mixed%nd = st_leafmtxp_scaled%nd
        st_leafmtxp_mixed%nlf = st_leafmtxp_scaled%nlf
        st_leafmtxp_mixed%nlfkt = st_leafmtxp_scaled%nlfkt
        st_leafmtxp_mixed%ktmax = st_leafmtxp_scaled%ktmax
        allocate(st_leafmtxp_mixed%st_lf(st_leafmtxp_mixed%nlf))

        ! 2. Check how many rows/columns we need for different data types, get the split values, then copy
        do ip = 1, st_leafmtxp_mixed%nlf
            SPLIT_CRITERION = 1e+1  ! CHANGE THIS HYPERPARAMETER!
            st_leafmtxp_mixed%st_lf(ip)%ltmtx = st_leafmtxp_scaled%st_lf(ip)%ltmtx  ! redundant I know
            st_leafmtxp_mixed%st_lf(ip)%ndl = st_leafmtxp_scaled%st_lf(ip)%ndl
            st_leafmtxp_mixed%st_lf(ip)%ndt = st_leafmtxp_scaled%st_lf(ip)%ndt
            st_leafmtxp_mixed%st_lf(ip)%nstrtl = st_leafmtxp_scaled%st_lf(ip)%nstrtl
            st_leafmtxp_mixed%st_lf(ip)%nstrtt = st_leafmtxp_scaled%st_lf(ip)%nstrtt

            if (st_leafmtxp_mixed%st_lf(ip)%ltmtx == 1) then
                st_leafmtxp_mixed%st_lf(ip)%kt = st_leafmtxp_scaled%st_lf(ip)%kt

                ! 1. need to find the max for relativity
                max = 0
                do ikt = 1, st_leafmtxp_mixed%st_lf(ip)%kt
                    if(st_leafmtxp_scaled%st_lf(ip)%diag(ikt) > max) then
                        max = st_leafmtxp_scaled%st_lf(ip)%diag(ikt)
                    endif
                enddo
                SPLIT_CRITERION = max * SPLIT_CRITERION! this becomes a reused SPLIT_CRITERION that accounts for relativity to max

                ! 2. now we know the relative SPLIT_CRITERION, we can count how many rows/columns should be 64-bit
                split_value = 0
                do ikt = 1, st_leafmtxp_mixed%st_lf(ip)%kt
                    if(st_leafmtxp_scaled%st_lf(ip)%diag(ikt) >= SPLIT_CRITERION) then
                        split_value = split_value + 1
                    endif
                enddo

                ! check the number of rows/cols that is suppose to be 64b
                ! print *, "split_value:"
                ! print *, ip, ": ", split_value

                st_leafmtxp_mixed%st_lf(ip)%kt_64 = split_value
                st_leafmtxp_mixed%st_lf(ip)%kt_32 = st_leafmtxp_mixed%st_lf(ip)%kt - split_value

                ! 3. Copy to new data type
                allocate(st_leafmtxp_mixed%st_lf(ip)%a1_64(st_leafmtxp_mixed%st_lf(ip)%ndt, st_leafmtxp_mixed%st_lf(ip)%kt_64))
                allocate(st_leafmtxp_mixed%st_lf(ip)%a1_32(st_leafmtxp_mixed%st_lf(ip)%ndt, st_leafmtxp_mixed%st_lf(ip)%kt_32))

                allocate(st_leafmtxp_mixed%st_lf(ip)%a2_64(st_leafmtxp_mixed%st_lf(ip)%ndl, st_leafmtxp_mixed%st_lf(ip)%kt_64))
                allocate(st_leafmtxp_mixed%st_lf(ip)%a2_32(st_leafmtxp_mixed%st_lf(ip)%ndl, st_leafmtxp_mixed%st_lf(ip)%kt_32))

                allocate(st_leafmtxp_mixed%st_lf(ip)%diag(st_leafmtxp_mixed%st_lf(ip)%kt))

                track_64 = 1
                track_32 = 1
                do ikt = 1, st_leafmtxp_mixed%st_lf(ip)%kt
                    if(st_leafmtxp_scaled%st_lf(ip)%diag(ikt) >= SPLIT_CRITERION) then
                        st_leafmtxp_mixed%st_lf(ip)%a1_64(1:st_leafmtxp_mixed%st_lf(ip)%ndt, track_64) = st_leafmtxp_scaled%st_lf(ip)%a1(1:st_leafmtxp_scaled%st_lf(ip)%ndt, ikt)
                        st_leafmtxp_mixed%st_lf(ip)%a2_64(1:st_leafmtxp_mixed%st_lf(ip)%ndl, track_64) = st_leafmtxp_scaled%st_lf(ip)%a2(1:st_leafmtxp_scaled%st_lf(ip)%ndl, ikt)
                        track_64 = track_64 + 1
                    else
                        st_leafmtxp_mixed%st_lf(ip)%a1_32(1:st_leafmtxp_mixed%st_lf(ip)%ndt, track_32) = st_leafmtxp_scaled%st_lf(ip)%a1(1:st_leafmtxp_scaled%st_lf(ip)%ndt, ikt)
                        st_leafmtxp_mixed%st_lf(ip)%a2_32(1:st_leafmtxp_mixed%st_lf(ip)%ndl, track_32) = st_leafmtxp_scaled%st_lf(ip)%a2(1:st_leafmtxp_scaled%st_lf(ip)%ndl, ikt)
                        track_32 = track_32 + 1
                    endif
                    st_leafmtxp_mixed%st_lf(ip)%diag(ikt) = st_leafmtxp_scaled%st_lf(ip)%diag(ikt)
                enddo

            elseif (st_leafmtxp_mixed%st_lf(ip)%ltmtx == 2) then
                ! TOASK: ARE FULL RANKS 100% 64-BITS???
                allocate(st_leafmtxp_mixed%st_lf(ip)%a1_64(st_leafmtxp_mixed%st_lf(ip)%ndt, st_leafmtxp_mixed%st_lf(ip)%ndl))
                st_leafmtxp_mixed%st_lf(ip)%a1_64(1:st_leafmtxp_mixed%st_lf(ip)%ndt, 1:st_leafmtxp_mixed%st_lf(ip)%ndl) = st_leafmtxp_scaled%st_lf(ip)%a1(1:st_leafmtxp_scaled%st_lf(ip)%ndt, 1:st_leafmtxp_scaled%st_lf(ip)%ndl)

            endif
        enddo

        ! ============================================== END OF PHASE 3 ================================================

! $omp critical
!         print *, "======= DEBUG ======="
!         do ip = 1, st_leafmtxp_32%nlf
!             print *, "YES"
!         enddo
!         print *, "mythread: ", omp_get_thread_num()
!         print *, st_leafmtxp%nd, " - ", st_leafmtxp_32%nd
!         print *, st_leafmtxp%nlf, " - ", st_leafmtxp_32%nlf
!         print *, st_leafmtxp%nlfkt, " - ", st_leafmtxp_32%nlfkt
!         print *, st_leafmtxp%ktmax, " - ", st_leafmtxp_32%ktmax
!         print *, "nd: ", nd
!         print *, "first five zu_32: ", zu_32(:5)

!         do ip = 1, 1
!             print *, st_leafmtxp%st_lf(ip)%ltmtx, " - ", st_leafmtxp_32%st_lf(ip)%ltmtx
!             ! print *, st_leafmtxp%st_lf(ip)%kt, " - ", st_leafmtxp_32%st_lf(ip)%kt
!             print *, st_leafmtxp%st_lf(ip)%ndl, " - ", st_leafmtxp_32%st_lf(ip)%ndl
!             print *, st_leafmtxp%st_lf(ip)%ndt, " - ", st_leafmtxp_32%st_lf(ip)%ndt
!             print *, st_leafmtxp%st_lf(ip)%nstrtl, " - ", st_leafmtxp_32%st_lf(ip)%nstrtl
!             print *, st_leafmtxp%st_lf(ip)%nstrtt, " - ", st_leafmtxp_32%st_lf(ip)%nstrtt
!         enddo
!         print *, "======= DEBUG ======="
! $omp end critical

! $omp single
        ! do ip = 1, 10
        !     print *, "== ip: ", ip, " =="
        !     print *, st_leafmtxp%st_lf(ip)%ltmtx, " - ", st_leafmtxp_32%st_lf(ip)%ltmtx
        !     print *, st_leafmtxp%st_lf(ip)%kt, " - ", st_leafmtxp_32%st_lf(ip)%kt
        !     print *, st_leafmtxp%st_lf(ip)%ndl, " - ", st_leafmtxp_32%st_lf(ip)%ndl
        !     print *, st_leafmtxp%st_lf(ip)%ndt, " - ", st_leafmtxp_32%st_lf(ip)%ndt
        !     print *, st_leafmtxp%st_lf(ip)%nstrtl, " - ", st_leafmtxp_32%st_lf(ip)%nstrtl
        !     print *, st_leafmtxp%st_lf(ip)%nstrtt, " - ", st_leafmtxp_32%st_lf(ip)%nstrtt
        ! enddo
! $omp end single
! $omp barrier

        lpmd => st_ctl%lpmd(:)
        lnp(0:) => st_ctl%lnp
        mpinr = lpmd(3)
        nrank = lpmd(2)
        icomm = lpmd(1)

        call MPI_Barrier(icomm, ierr)
        ! st_measure_time = MPI_Wtime()

        mstep = param(83)
        eps = param(91)
        allocate(wws(maxval(lnp(0:nrank-1))), wwr(maxval(lnp(0:nrank-1))))
        allocate(zr(nd), zshdw(nd), zp(nd), zt(nd), zkp(nd), zakp(nd), zkt(nd), zakt(nd))
        bnorm = dsqrt(HACApK_dotp_d(nd, b, b))

        if(st_ctl%param(1)>0 .and. mpinr==0)    print *, 'HACApK_bicgstab_lfmtx_hyp start'

!$omp parallel

!$omp workshare
        zp(1:nd) = 0.0d0
        zakp(1:nd) = 0.0d0
        zr(:nd) = b(:nd)
!$omp end workshare

        ! =============== matvecmul timer starts here =================
        ! Everything here has no relationship to bicgstab
        ! COMMENT THEM ALL OUT WHEN WANT TO USE BICGSTAB

        ! === Rigorous Testing: 1000 inner loops sandwiched between omp master ===
        ! call HACApK_adot_lfmtx_hyp(zshdw, st_leafmtxp, st_ctl, dummy_64, wws, wwr, isct, irct, nd)
        ! call HACApK_adot_lfmtx_hyp_32_64(zshdw, st_leafmtxp_32, st_ctl, dummy_64, wws, wwr, isct, irct, nd)
        ! call HACApK_adot_lfmtx_hyp_32(zshdw, st_leafmtxp_32, st_ctl, dummy_32, wws, wwr, isct, irct, nd)

        ! === Scaled Leaves' Rigorous Testing: 1000 inner loops sandwiched between omp master ===
        ! call HACApK_adot_lfmtx_hyp_scaled(zshdw, st_leafmtxp_scaled, st_ctl, dummy_64, wws, wwr, isct, irct, nd)
        ! call HACApK_adot_lfmtx_hyp_scaled_32_64(zshdw, st_leafmtxp_scaled_32, st_ctl, dummy_64, wws, wwr, isct, irct, nd)
        ! call HACApK_adot_lfmtx_hyp_scaled_32(zshdw, st_leafmtxp_scaled_32, st_ctl, dummy_32, wws, wwr, isct, irct, nd)

        ! === Mixed Precision Rigorous Testing with 1000 inner loops
        call HACApK_adot_lfmtx_hyp_mixed(zshdw, st_leafmtxp_mixed, st_ctl, dummy_64, wws, wwr, isct, irct, nd)
        ! =============== matvecmul timer stops here =================
        stop

        ! everything below is usual bicgstab
        ! There are three functions to call, one here, the other two calculating zakp and zakt, totaling 3
        ! Also means 3 placeholder_32

        ! I need more information on this function
        ! === Phase 1
        ! call HACApK_adotsub_lfmtx_hyp(zr, zshdw, st_leafmtxp, st_ctl, u, wws, wwr, isct, irct, nd)
        ! call HACApK_adotsub_lfmtx_hyp_32_64(zr, zshdw, st_leafmtxp_32, st_ctl, u, wws, wwr, isct, irct, nd)
        ! placeholder_32(:nd) = u(:nd)
        ! call HACApK_adotsub_lfmtx_hyp_32(zr, zshdw, st_leafmtxp_32, st_ctl, placeholder_32, wws, wwr, isct, irct, nd)

        ! === Phase 2: Scaled Leaves
        ! call HACApK_adotsub_lfmtx_hyp_scaled(zr, zshdw, st_leafmtxp_scaled, st_ctl, u, wws, wwr, isct, irct, nd)
        ! call HACApK_adotsub_lfmtx_hyp_scaled_32_64(zr, zshdw, st_leafmtxp_scaled_32, st_ctl, u, wws, wwr, isct, irct, nd)
        ! placeholder_32(:nd) = u(:nd)
        ! call HACApK_adotsub_lfmtx_hyp_scaled_32(zr, zshdw, st_leafmtxp_scaled_32, st_ctl, placeholder_32, wws, wwr, isct, irct, nd)

        ! === Phase 3: Mixed Data Structure
        ! call HACApK_adotsub_lfmtx_hyp_mixed(zr, zshdw, st_leafmtxp_mixed, st_ctl, u, wws, wwr, isct, irct, nd)

!$omp barrier

!$omp workshare
        zshdw(:nd) = zr(:nd)
!$omp end workshare

!$omp single
        zrnorm = dsqrt(HACApK_dotp_d(nd, zr, zr))
        if(mpinr==0)        print *, 'Original relative residual norm =', zrnorm/bnorm, bnorm
!$omp end single

        !---------------------------------------------------------
        !---------------- is here the loop source? ---------------
        !---------------------------------------------------------
        !----- mstep(param(83)) is initialized in HACApK_init, ---
        !-- renewed in  "st_ctl%param(83) = ppohBEM_max_steps" ---
        !------------- in the bem_bb_fw_HACApK_0.4.1. ------------
        !---------------------------------------------------------
        !--- ppohBEM_max_steps is set in Read_bem_bb_config, -----
        !------------ written in bem-bb-config.txt. --------------
        !---------------------------------------------------------

! === bicgstab timing starts here ===
!$omp barrier
!$omp master
        st_measure_time = MPI_Wtime()
!$omp end master
! === ===

        mstep = 150  ! manual maximum steps allowed
        do in = 1, mstep
        ! do in = 1, 2
            if(zrnorm/bnorm < eps) then
!            if(zrnorm/bnorm < 100d0) then
!$omp barrier
!$omp master
              en_measure_time = MPI_Wtime()
              print *, "Time taken for bicgstab to converge:"
              print *, en_measure_time - st_measure_time
              print *, "End of Bicgstab"
!$omp end master

    ! ! DEBUG U
    ! !$omp master
    ! print *, '**************** DEBUG U ****************'

    ! open(10, file='debug_u.data')
    !   do i = 1, nd
    !     write(10, *) u(i)
    !   enddo
    ! close(10)

    ! !$omp end master
    ! !$omp barrier

              exit
            endif

!$omp workshare
            zp(:nd)  = zr(:nd) + beta * (zp(:nd)-zeta * zakp(:nd))
            zkp(:nd) = zp(:nd)
            placeholder_32(:nd) = zkp(:nd)  ! for 32 bit
!$omp end workshare

            ! call HACApK_adot_lfmtx_hyp(zakp, st_leafmtxp, st_ctl, zkp, wws, wwr, isct, irct, nd)
            ! call HACApK_adot_lfmtx_hyp_32_64(zakp, st_leafmtxp_32, st_ctl, zkp, wws, wwr, isct, irct, nd)
            ! remember to set placeholder before calling pure 32!
            ! call HACApK_adot_lfmtx_hyp_32(zakp, st_leafmtxp_32, st_ctl, placeholder_32, wws, wwr, isct, irct, nd)

            ! call HACApK_adot_lfmtx_hyp_scaled(zakp, st_leafmtxp_scaled, st_ctl, zkp, wws, wwr, isct, irct, nd)
            ! call HACApK_adot_lfmtx_hyp_scaled_32_64(zakp, st_leafmtxp_scaled_32, st_ctl, zkp, wws, wwr, isct, irct, nd)
            ! remember to set placeholder before calling pure 32!
            ! call HACApK_adot_lfmtx_hyp_scaled_32(zakp, st_leafmtxp_scaled_32, st_ctl, placeholder_32, wws, wwr, isct, irct, nd)

            ! call HACApK_adot_lfmtx_hyp_mixed(zakp, st_leafmtxp_mixed, st_ctl, zkp, wws, wwr, isct, irct, nd)

!$omp barrier

!$omp single
            znom  = HACApK_dotp_d(nd, zshdw, zr)
            alpha = znom / HACApK_dotp_d(nd, zshdw, zakp)       !znom/zden
            znomold = znom
!$omp end single

!$omp workshare
            zt(:nd)  = zr(:nd) - alpha * zakp(:nd)
            zkt(:nd) = zt(:nd)
            placeholder_32(:nd) = zkt(:nd)  ! for 32 bit
!$omp end workshare

            ! call HACApK_adot_lfmtx_hyp(zakt, st_leafmtxp, st_ctl, zkt, wws, wwr, isct, irct, nd)
            ! call HACApK_adot_lfmtx_hyp_32_64(zakt, st_leafmtxp_32, st_ctl, zkt, wws, wwr, isct, irct, nd)
            ! remember to set placeholder before calling pure 32!
            ! call HACApK_adot_lfmtx_hyp_32(zakt, st_leafmtxp_32, st_ctl, placeholder_32, wws, wwr, isct, irct, nd)

            ! call HACApK_adot_lfmtx_hyp_scaled(zakt, st_leafmtxp_scaled, st_ctl, zkt, wws, wwr, isct, irct, nd)
            ! call HACApK_adot_lfmtx_hyp_scaled_32_64(zakt, st_leafmtxp_scaled_32, st_ctl, zkt, wws, wwr, isct, irct, nd)
            ! remember to set placeholder before calling pure 32!
            ! call HACApK_adot_lfmtx_hyp_scaled_32(zakt, st_leafmtxp_scaled_32, st_ctl, placeholder_32, wws, wwr, isct, irct, nd)

            ! call HACApK_adot_lfmtx_hyp_mixed(zakt, st_leafmtxp_mixed, st_ctl, zkt, wws, wwr, isct, irct, nd)

!$omp barrier

!$omp single
            znom = HACApK_dotp_d(nd, zakt, zt)
            zeta = znom / HACApK_dotp_d(nd, zakt, zakt)         ! znom/zden
!$omp end single

!$omp workshare
            u(:nd)  = u(:nd) + alpha * zkp(:nd) + zeta * zkt(:nd)  ! why is this u updated if never been used?
            zr(:nd) = zt(:nd) - zeta * zakt(:nd)
!$omp end workshare

!$omp single
            beta = alpha/zeta * HACApK_dotp_d(nd, zshdw, zr) / znomold
            zrnorm = dsqrt(HACApK_dotp_d(nd, zr, zr))
            nstp = in
            call MPI_Barrier(icomm, ierr)

            ! en_measure_time = MPI_Wtime()
            ! time = en_measure_time - st_measure_time

            ! if(st_ctl%param(1)>0 .and. mpinr==0)        print *, in, time, log10(zrnorm/bnorm)
            if(st_ctl%param(1)>0 .and. mpinr==0)        print *, in, log10(zrnorm/bnorm)
!$omp end single
        enddo

        ! time_total = time_total + time
!       if(getIDinAllProcs()==1)        print *, 'total time', time_total

! Added by Iwashita & Rise on 30 November 2018
! To check the recalculated relative resideual norm as compared to the original.
! Weird bug

!$omp workshare
!zr(:nd)=b(:nd)
!zshdw(:nd)=0d0
  zakp(:nd)=0d0
  zr(:nd)=0d0
!$omp end workshare

!call HACApK_adotsub_lfmtx_hyp(zr,zshdw,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
 call HACApK_adot_lfmtx_hyp(zakp,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
!$omp barrier

!$omp single
  do i=1,nd
    zr(i)=b(i)-zakp(i)
  enddo

  zrnorm=HACApK_dotp_d(nd,zr,zr)
  zrnorm=dsqrt(zrnorm)
!write(6,*) dsqrt(zrnorm)/bnorm
    print *, "Recalculated relative residual norm = ", zrnorm/bnorm, bnorm

!$omp end single

!!$omp end parallel
!     open (10,file="solx.data")
!     do i=1,nd
!       write(10,*) u(i)
!     enddo
!     close(10)


!!$omp single
!        zr(1:nd)=b(1:nd)
!  !$omp end single
!   !$omp barrier
!        call HACApK_adotsub_lfmtx_hyp(zr, zshdw, st_leafmtxp, st_ctl, u, wws, wwr, isct, irct, nd)
!!$omp barrier
!!$omp single
!        zrnorm = dsqrt(HACApK_dotp_d(nd, zr, zr))
!        print *, "Recalculated relative residual norm = ", zrnorm/bnorm, bnorm
!!$omp end single


!$omp end parallel

    ! Added by Rise
    ! print *, 'total time', time_total

    end subroutine HACApK_bicgstab_lfmtx_hyp


!***HACApK_adotsub_lfmtx_hyp (Original Version)
    subroutine HACApK_adotsub_lfmtx_hyp(zr, zau, st_leafmtxp, st_ctl, zu, wws, wwr, isct, irct, nd)
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zr(*), zau(*), zu(*), wws(*), wwr(*)
        integer*4 :: isct(*), irct(*)

        ! print *, zu(:5)
        ! stop "debug"

        call HACApK_adot_lfmtx_hyp(zau, st_leafmtxp, st_ctl, zu, wws, wwr, isct, irct, nd)

!$omp barrier
!$omp workshare
        zr(1:nd) = zr(1:nd)-zau(1:nd)
!$omp end workshare
    end subroutine HACApK_adotsub_lfmtx_hyp


!***HACApK_adotsub_lfmtx_hyp (Rise's version)
    subroutine HACApK_adotsub_lfmtx_hyp_32(zr, zau, st_leafmtxp_32, st_ctl, zu, wws, wwr, isct, irct, nd)
        type(st_HACApK_leafmtxp_32) :: st_leafmtxp_32
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zr(*), zau(*), wws(*), wwr(*)
        real*4 :: zu(*)
        integer*4 :: isct(*), irct(*)

        call HACApK_adot_lfmtx_hyp_32(zau, st_leafmtxp_32, st_ctl, zu, wws, wwr, isct, irct, nd)


!$omp barrier
!$omp workshare
        zr(1:nd) = zr(1:nd)-zau(1:nd)
!$omp end workshare
    end subroutine HACApK_adotsub_lfmtx_hyp_32


!***HACApK_adotsub_lfmtx_hyp (Rise's version)
    subroutine HACApK_adotsub_lfmtx_hyp_32_64(zr, zau, st_leafmtxp_32, st_ctl, zu, wws, wwr, isct, irct, nd)
        type(st_HACApK_leafmtxp_32) :: st_leafmtxp_32
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zr(*), zau(*), wws(*), wwr(*)
        real*8 :: zu(*)
        integer*4 :: isct(*), irct(*)

        call HACApK_adot_lfmtx_hyp_32_64(zau, st_leafmtxp_32, st_ctl, zu, wws, wwr, isct, irct, nd)


!$omp barrier
!$omp workshare
        zr(1:nd) = zr(1:nd)-zau(1:nd)
!$omp end workshare
    end subroutine HACApK_adotsub_lfmtx_hyp_32_64


! ============ SCALED VERSIONS =============
!***HACApK_adotsub_lfmtx_hyp (Original Version)
    subroutine HACApK_adotsub_lfmtx_hyp_scaled(zr, zau, st_leafmtxp_scaled, st_ctl, zu, wws, wwr, isct, irct, nd)
        type(st_HACApK_leafmtxp_scaled) :: st_leafmtxp_scaled
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zr(*), zau(*), zu(*), wws(*), wwr(*)
        integer*4 :: isct(*), irct(*)

        ! print *, zu(:5)
        ! stop "debug"

        call HACApK_adot_lfmtx_hyp_scaled(zau, st_leafmtxp_scaled, st_ctl, zu, wws, wwr, isct, irct, nd)

!$omp barrier
!$omp workshare
        zr(1:nd) = zr(1:nd)-zau(1:nd)
!$omp end workshare
    end subroutine HACApK_adotsub_lfmtx_hyp_scaled


!***HACApK_adotsub_lfmtx_hyp (Rise's version)
    subroutine HACApK_adotsub_lfmtx_hyp_scaled_32(zr, zau, st_leafmtxp_scaled_32, st_ctl, zu, wws, wwr, isct, irct, nd)
        type(st_HACApK_leafmtxp_scaled_32) :: st_leafmtxp_scaled_32
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zr(*), zau(*), wws(*), wwr(*)
        real*4 :: zu(*)
        integer*4 :: isct(*), irct(*)

        call HACApK_adot_lfmtx_hyp_scaled_32(zau, st_leafmtxp_scaled_32, st_ctl, zu, wws, wwr, isct, irct, nd)


!$omp barrier
!$omp workshare
        zr(1:nd) = zr(1:nd)-zau(1:nd)
!$omp end workshare
    end subroutine HACApK_adotsub_lfmtx_hyp_scaled_32


!***HACApK_adotsub_lfmtx_hyp (Rise's version)
    subroutine HACApK_adotsub_lfmtx_hyp_scaled_32_64(zr, zau, st_leafmtxp_scaled_32, st_ctl, zu, wws, wwr, isct, irct, nd)
        type(st_HACApK_leafmtxp_scaled_32) :: st_leafmtxp_scaled_32
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zr(*), zau(*), wws(*), wwr(*)
        real*8 :: zu(*)
        integer*4 :: isct(*), irct(*)

        call HACApK_adot_lfmtx_hyp_scaled_32_64(zau, st_leafmtxp_scaled_32, st_ctl, zu, wws, wwr, isct, irct, nd)


!$omp barrier
!$omp workshare
        zr(1:nd) = zr(1:nd)-zau(1:nd)
!$omp end workshare
    end subroutine HACApK_adotsub_lfmtx_hyp_scaled_32_64


subroutine HACApK_adotsub_lfmtx_hyp_mixed(zr, zau, st_leafmtxp_mixed, st_ctl, zu, wws, wwr, isct, irct, nd)
        type(st_HACApK_leafmtxp_mixed) :: st_leafmtxp_mixed
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zr(*), zau(*), zu(*), wws(*), wwr(*)
        integer*4 :: isct(*), irct(*)

        ! print *, zu(:5)
        ! stop "debug"

        call HACApK_adot_lfmtx_hyp_mixed(zau, st_leafmtxp_mixed, st_ctl, zu, wws, wwr, isct, irct, nd)

!$omp barrier
!$omp workshare
        zr(1:nd) = zr(1:nd)-zau(1:nd)
!$omp end workshare
    end subroutine HACApK_adotsub_lfmtx_hyp_mixed


!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!***HACApK_gcrm_lfmtx
    subroutine HACApK_gcrm_lfmtx(st_leafmtxp, st_ctl, st_bemv, u, b, param, nd, nstp, lrtrn)
        include 'mpif.h'
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl
        type(st_HACApK_calc_entry) :: st_bemv
        real*8 :: u(nd), b(nd)
        real*8 :: param(*)

        real*8, dimension(:), allocatable :: zr, zar, capap
        real*8, dimension(:, :), allocatable, target :: zp, zap
        real*8, pointer :: zq(:)
        real*8, dimension(:), allocatable :: wws, wwr
        integer*4 :: isct(2), irct(2)
        integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)

        lpmd => st_ctl%lpmd(:)
        lnp(0:) => st_ctl%lnp
        mpinr = lpmd(3)
        nrank = lpmd(2)
        icomm = lpmd(1)

        call MPI_Barrier(icomm, ierr)
        st_measure_time = MPI_Wtime()
        if(st_ctl%param(1)>0 .and. mpinr==0)        print *, 'gcr_lfmtx_hyp start'

        mstep = param(83)
        mreset = param(87)
        eps = param(91)
        allocate(wws(maxval(lnp(0:nrank-1))), wwr(maxval(lnp(0:nrank-1))))
        allocate(zr(nd), zar(nd), zp(nd, mreset), zap(nd, mreset), capap(mreset))
        alpha = 0.0
        bnorm = dsqrt(HACApK_dotp_d(nd, b, b))
        call HACApK_adot_lfmtx_hyp(zar, st_leafmtxp, st_ctl, u, wws, wwr, isct, irct, nd)

        zr(:nd) = b(:nd) - zar(:nd)
        zp(:nd, 1) = zr(:nd)
        zrnorm = dsqrt(HACApK_dotp_d(nd, zr, zr))

        call MPI_Barrier( icomm, ierr )
        en_measure_time = MPI_Wtime()
        time = en_measure_time - st_measure_time
        if(st_ctl%param(1)>0 .and. mpinr==0)        print *, 0, time, log10(zrnorm/bnorm)
        if(zrnorm/bnorm < eps)      return

        call HACApK_adot_lfmtx_hyp(zap(:nd, 1), st_leafmtxp, st_ctl, zp(:nd, 1), wws, wwr, isct, irct, nd)

        do in = 1, mstep
            ik = mod(in-1, mreset) + 1
            zq => zap(:nd, ik)

            znom = HACApK_dotp_d(nd, zq, zr)
            capap(ik) = HACApK_dotp_d(nd, zq, zq)
            alpha = znom/capap(ik)
            u(:nd)  = u(:nd)  + alpha * zp(:nd, ik)
            zr(:nd) = zr(:nd) - alpha * zq(:nd)
            zrnorm = dsqrt(HACApK_dotp_d(nd, zr, zr))

            call MPI_Barrier(icomm, ierr)
            en_measure_time = MPI_Wtime()
            time = en_measure_time - st_measure_time
            if(st_ctl%param(1)>0 .and. mpinr==0)        print *, in, time, log10(zrnorm/bnorm)

            if(zrnorm/bnorm < eps .or. in==mstep)       exit

            call HACApK_adot_lfmtx_hyp(zar, st_leafmtxp, st_ctl, zr, wws, wwr, isct, irct, nd)
            ikn = mod(in, mreset) + 1
            zp(:nd, ikn)  = zr(:nd)
            zap(:nd, ikn) = zar(:nd)

            do il = 1, ik
                zq => zap(:nd, il)
                znom = HACApK_dotp_d(nd, zq, zar)
                beta = -znom/capap(il)
                zp(:nd, ikn)  = zp(:nd, ikn)  + beta * zp(:nd, il)
                zap(:nd, ikn) = zap(:nd, ikn) + beta * zq(:nd)
            enddo
        enddo
        nstp = in
    end subroutine


 !***HACApK_adot_lfmtx_hyp_mixed (PHASE 3: MIXED DATA STRUCTURE VERSION)
subroutine HACApK_adot_lfmtx_hyp_mixed(zau, st_leafmtxp_mixed, st_ctl, zu, wws, wwr, isct, irct, nd)
    include 'mpif.h'
    type(st_HACApK_leafmtxp_mixed) :: st_leafmtxp_mixed
    type(st_HACApK_lcontrol) :: st_ctl
    real*8 :: zau(*), wws(*), wwr(*)
    real*8 :: zu(*)
    integer*4 :: isct(*), irct(*)
    integer*4 :: ISTATUS(MPI_STATUS_SIZE)
    integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
    real*8 :: time_start, time_end = 0.0
    ! integer*4, dimension(:), allocatable :: ISTATUS
    1000 format(5(a, i10)/)
    2000 format(5(a, f10.4)/)

    lpmd => st_ctl%lpmd(:)
    lnp(0:) => st_ctl%lnp
    lsp(0:) => st_ctl%lsp
    lthr(0:) => st_ctl%lthr
    ! allocate(ISTATUS(MPI_STATUS_SIZE))
    mpinr = lpmd(3)
    mpilog = lpmd(4)
    nrank = lpmd(2)
    icomm = lpmd(1)
    ndnr_s = lpmd(6)
    ndnr_e = lpmd(7)
    ndnr = lpmd(5)
    zau(:nd) = 0.0d0

    !$omp barrier

    !$omp master
      time_start = MPI_Wtime()
    !$omp end master

    do i = 1, 1000
      call HACApK_adot_body_lfmtx_hyp_mixed(zau, st_leafmtxp_mixed, st_ctl, zu, nd)
    enddo

    !$omp barrier

    !$omp master
      time_end = MPI_Wtime()
      print *, "Total execution time:"
      print *, time_end - time_start
    !$omp end master

    !$omp barrier

    !=== printing part
    ! !$omp master
    ! ! print *, 'DEBUG 64'

    ! open(10, file='result_vector_mixed_1e0.data')
    !   do i = 1, nd
    !     write(10, *) zau(i)
    !   enddo
    ! close(10)

    ! !$omp end master
    ! !$omp barrier

    stop

    !$omp master
    if (nrank > 1) then
        wws(1:lnp(mpinr)) = zau(lsp(mpinr):lsp(mpinr) + lnp(mpinr) - 1)
        ncdp = mod(mpinr + 1, nrank)
        ncsp = mod(mpinr + nrank - 1, nrank)
        isct(1) = lnp(mpinr);isct(2) = lsp(mpinr)
        do ic = 1,nrank-1
            call MPI_SENDRECV(isct, 2, MPI_INTEGER, ncdp, 1, irct,2, MPI_INTEGER, ncsp, 1, icomm, ISTATUS, ierr)
            call MPI_SENDRECV(wws, isct, MPI_DOUBLE_PRECISION, ncdp, 1, wwr, irct, MPI_DOUBLE_PRECISION, ncsp, 1, icomm, ISTATUS, ierr)
            zau(irct(2):irct(2) + irct(1) - 1) = zau(irct(2):irct(2) + irct(1) - 1) + wwr(:irct(1))
            wws(:irct(1)) = wwr(:irct(1))
            isct(:2) = irct(:2)
        enddo
    endif
    !$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp_mixed



 !***HACApK_adot_lfmtx_hyp_scaled (Scaled version)
subroutine HACApK_adot_lfmtx_hyp_scaled(zau, st_leafmtxp_scaled, st_ctl, zu, wws, wwr, isct, irct, nd)
    include 'mpif.h'
    type(st_HACApK_leafmtxp_scaled) :: st_leafmtxp_scaled
    type(st_HACApK_lcontrol) :: st_ctl
    real*8 :: zau(*), wws(*), wwr(*)
    real*8 :: zu(*)
    integer*4 :: isct(*), irct(*)
    integer*4 :: ISTATUS(MPI_STATUS_SIZE)
    integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
    real*8 :: time_start, time_end = 0.0
    ! integer*4, dimension(:), allocatable :: ISTATUS
    1000 format(5(a, i10)/)
    2000 format(5(a, f10.4)/)

    lpmd => st_ctl%lpmd(:)
    lnp(0:) => st_ctl%lnp
    lsp(0:) => st_ctl%lsp
    lthr(0:) => st_ctl%lthr
    ! allocate(ISTATUS(MPI_STATUS_SIZE))
    mpinr = lpmd(3)
    mpilog = lpmd(4)
    nrank = lpmd(2)
    icomm = lpmd(1)
    ndnr_s = lpmd(6)
    ndnr_e = lpmd(7)
    ndnr = lpmd(5)
    zau(:nd) = 0.0d0

    ! !$omp barrier

    ! !$omp master
    !   time_start = MPI_Wtime()
    ! !$omp end master

    ! do i = 1, 1000
      call HACApK_adot_body_lfmtx_hyp_scaled_64(zau, st_leafmtxp_scaled, st_ctl, zu, nd)
    ! enddo

    ! !$omp barrier

    ! !$omp master
    !   time_end = MPI_Wtime()
    !   print *, "Total execution time:"
    !   print *, time_end - time_start
    ! !$omp end master

    ! !$omp barrier

    ! !$omp master
    ! ! print *, 'DEBUG 64'

    ! open(10, file='result_vector_scaled_64.data')
    !   do i = 1, nd
    !     write(10, *) zau(i)
    !   enddo
    ! close(10)

    ! !$omp end master
    ! !$omp barrier

    ! stop

    !$omp master
    if (nrank > 1) then
        wws(1:lnp(mpinr)) = zau(lsp(mpinr):lsp(mpinr) + lnp(mpinr) - 1)
        ncdp = mod(mpinr + 1, nrank)
        ncsp = mod(mpinr + nrank - 1, nrank)
        isct(1) = lnp(mpinr);isct(2) = lsp(mpinr)
        do ic = 1,nrank-1
            call MPI_SENDRECV(isct, 2, MPI_INTEGER, ncdp, 1, irct,2, MPI_INTEGER, ncsp, 1, icomm, ISTATUS, ierr)
            call MPI_SENDRECV(wws, isct, MPI_DOUBLE_PRECISION, ncdp, 1, wwr, irct, MPI_DOUBLE_PRECISION, ncsp, 1, icomm, ISTATUS, ierr)
            zau(irct(2):irct(2) + irct(1) - 1) = zau(irct(2):irct(2) + irct(1) - 1) + wwr(:irct(1))
            wws(:irct(1)) = wwr(:irct(1))
            isct(:2) = irct(:2)
        enddo
    endif
    !$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp_scaled


 subroutine HACApK_adot_lfmtx_hyp_scaled_32(zau, st_leafmtxp_scaled_32, st_ctl, zu_32, wws, wwr, isct, irct, nd)
    include 'mpif.h'
    type(st_HACApK_leafmtxp_scaled_32) :: st_leafmtxp_scaled_32
    type(st_HACApK_lcontrol) :: st_ctl
    real*8 :: zau(*), wws(*), wwr(*)
    real*4 :: zu_32(*)
    integer*4 :: isct(*), irct(*)
    integer*4 :: ISTATUS(MPI_STATUS_SIZE)
    integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
    real*8 :: time_start, time_end = 0.0
    ! integer*4, dimension(:), allocatable :: ISTATUS
    1000 format(5(a, i10)/)
    2000 format(5(a, f10.4)/)

    lpmd => st_ctl%lpmd(:)
    lnp(0:) => st_ctl%lnp
    lsp(0:) => st_ctl%lsp
    lthr(0:) => st_ctl%lthr
    ! allocate(ISTATUS(MPI_STATUS_SIZE))
    mpinr = lpmd(3)
    mpilog = lpmd(4)
    nrank = lpmd(2)
    icomm = lpmd(1)
    ndnr_s = lpmd(6)
    ndnr_e = lpmd(7)
    ndnr = lpmd(5)
    zau(:nd) = 0.0d0

    ! !$omp barrier

    ! !$omp master
    !   time_start = MPI_Wtime()
    ! !$omp end master

    ! do i = 1, 1000
      call HACApK_adot_body_lfmtx_hyp_scaled_32(zau, st_leafmtxp_scaled_32, st_ctl, zu_32, nd)
    ! enddo

    ! !$omp barrier

    ! !$omp master
    !   time_end = MPI_Wtime()
    !   print *, "Total execution time:"
    !   print *, time_end - time_start
    ! !$omp end master

    ! !$omp barrier

    ! !$omp master
    ! print *, 'DEBUG 64'

    ! open(10, file='result_vector_scaled_32o32o32x32+64.data')
    !   do i = 1, nd
    !     write(10, *) zau(i)
    !   enddo
    ! close(10)

    ! !$omp end master
    ! !$omp barrier

    ! stop

    !$omp master
    if (nrank > 1) then
        wws(1:lnp(mpinr)) = zau(lsp(mpinr):lsp(mpinr) + lnp(mpinr) - 1)
        ncdp = mod(mpinr + 1, nrank)
        ncsp = mod(mpinr + nrank - 1, nrank)
        isct(1) = lnp(mpinr);isct(2) = lsp(mpinr)
        do ic = 1,nrank-1
            call MPI_SENDRECV(isct, 2, MPI_INTEGER, ncdp, 1, irct,2, MPI_INTEGER, ncsp, 1, icomm, ISTATUS, ierr)
            call MPI_SENDRECV(wws, isct, MPI_DOUBLE_PRECISION, ncdp, 1, wwr, irct, MPI_DOUBLE_PRECISION, ncsp, 1, icomm, ISTATUS, ierr)
            zau(irct(2):irct(2) + irct(1) - 1) = zau(irct(2):irct(2) + irct(1) - 1) + wwr(:irct(1))
            wws(:irct(1)) = wwr(:irct(1))
            isct(:2) = irct(:2)
        enddo
    endif
    !$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp_scaled_32


  subroutine HACApK_adot_lfmtx_hyp_scaled_32_64(zau, st_leafmtxp_scaled_32, st_ctl, zu, wws, wwr, isct, irct, nd)
    include 'mpif.h'
    type(st_HACApK_leafmtxp_scaled_32) :: st_leafmtxp_scaled_32
    type(st_HACApK_lcontrol) :: st_ctl
    real*8 :: zau(*), wws(*), wwr(*)
    real*8 :: zu(*)
    integer*4 :: isct(*), irct(*)
    integer*4 :: ISTATUS(MPI_STATUS_SIZE)
    integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
    real*8 :: time_start, time_end = 0.0
    ! integer*4, dimension(:), allocatable :: ISTATUS
    1000 format(5(a, i10)/)
    2000 format(5(a, f10.4)/)

    lpmd => st_ctl%lpmd(:)
    lnp(0:) => st_ctl%lnp
    lsp(0:) => st_ctl%lsp
    lthr(0:) => st_ctl%lthr
    ! allocate(ISTATUS(MPI_STATUS_SIZE))
    mpinr = lpmd(3)
    mpilog = lpmd(4)
    nrank = lpmd(2)
    icomm = lpmd(1)
    ndnr_s = lpmd(6)
    ndnr_e = lpmd(7)
    ndnr = lpmd(5)
    zau(:nd) = 0.0d0

    ! !$omp barrier

    ! !$omp master
    !   time_start = MPI_Wtime()
    ! !$omp end master

    ! do i = 1, 1000
      call HACApK_adot_body_lfmtx_hyp_scaled_mixed(zau, st_leafmtxp_scaled_32, st_ctl, zu, nd)
    ! enddo

    ! !$omp barrier

    ! !$omp master
    !   time_end = MPI_Wtime()
    !   print *, "Total execution time:"
    !   print *, time_end - time_start
    ! !$omp end master

    ! !$omp barrier

    ! !$omp master
    ! print *, 'DEBUG 64'

    ! open(10, file='result_vector_scaled_32o32o32x64+64.data')
    !   do i = 1, nd
    !     write(10, *) zau(i)
    !   enddo
    ! close(10)

    ! !$omp end master
    ! !$omp barrier

    ! stop

    !$omp master
    if (nrank > 1) then
        wws(1:lnp(mpinr)) = zau(lsp(mpinr):lsp(mpinr) + lnp(mpinr) - 1)
        ncdp = mod(mpinr + 1, nrank)
        ncsp = mod(mpinr + nrank - 1, nrank)
        isct(1) = lnp(mpinr);isct(2) = lsp(mpinr)
        do ic = 1,nrank-1
            call MPI_SENDRECV(isct, 2, MPI_INTEGER, ncdp, 1, irct,2, MPI_INTEGER, ncsp, 1, icomm, ISTATUS, ierr)
            call MPI_SENDRECV(wws, isct, MPI_DOUBLE_PRECISION, ncdp, 1, wwr, irct, MPI_DOUBLE_PRECISION, ncsp, 1, icomm, ISTATUS, ierr)
            zau(irct(2):irct(2) + irct(1) - 1) = zau(irct(2):irct(2) + irct(1) - 1) + wwr(:irct(1))
            wws(:irct(1)) = wwr(:irct(1))
            isct(:2) = irct(:2)
        enddo
    endif
    !$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp_scaled_32_64


 !***HACApK_adot_lfmtx_hyp (Original 64-bit version)
subroutine HACApK_adot_lfmtx_hyp(zau, st_leafmtxp, st_ctl, zu, wws, wwr, isct, irct, nd)
    include 'mpif.h'
    type(st_HACApK_leafmtxp) :: st_leafmtxp
    type(st_HACApK_lcontrol) :: st_ctl
    real*8 :: zau(*), wws(*), wwr(*)
    real*8 :: zu(*)
    integer*4 :: isct(*), irct(*)
    integer*4 :: ISTATUS(MPI_STATUS_SIZE)
    integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
    real*8 :: time_start, time_end = 0.0
    ! integer*4, dimension(:), allocatable :: ISTATUS
    1000 format(5(a, i10)/)
    2000 format(5(a, f10.4)/)

    lpmd => st_ctl%lpmd(:)
    lnp(0:) => st_ctl%lnp
    lsp(0:) => st_ctl%lsp
    lthr(0:) => st_ctl%lthr
    ! allocate(ISTATUS(MPI_STATUS_SIZE))
    mpinr = lpmd(3)
    mpilog = lpmd(4)
    nrank = lpmd(2)
    icomm = lpmd(1)
    ndnr_s = lpmd(6)
    ndnr_e = lpmd(7)
    ndnr = lpmd(5)
    zau(:nd) = 0.0d0

    ! !$omp barrier

    ! !$omp master
    !   time_start = MPI_Wtime()
    ! !$omp end master

    ! do i = 1, 1000
      call HACApK_adot_body_lfmtx_hyp_rise_64(zau, st_leafmtxp, st_ctl, zu, nd)
    ! enddo

    ! !$omp barrier

    ! !$omp master
    !   time_end = MPI_Wtime()
    !   print *, "Total execution time:"
    !   print *, time_end - time_start
    ! !$omp end master

    ! !$omp barrier

    ! !$omp master
    ! ! ! print *, 'DEBUG 64'

    ! open(10, file='result_vector_64.data')
    !   do i = 1, nd
    !     write(10, *) zau(i)
    !   enddo
    ! close(10)

    ! !$omp end master
    ! !$omp barrier

    ! stop

    !$omp master
    if (nrank > 1) then
        wws(1:lnp(mpinr)) = zau(lsp(mpinr):lsp(mpinr) + lnp(mpinr) - 1)
        ncdp = mod(mpinr + 1, nrank)
        ncsp = mod(mpinr + nrank - 1, nrank)
        isct(1) = lnp(mpinr);isct(2) = lsp(mpinr)
        do ic = 1,nrank-1
            call MPI_SENDRECV(isct, 2, MPI_INTEGER, ncdp, 1, irct,2, MPI_INTEGER, ncsp, 1, icomm, ISTATUS, ierr)
            call MPI_SENDRECV(wws, isct, MPI_DOUBLE_PRECISION, ncdp, 1, wwr, irct, MPI_DOUBLE_PRECISION, ncsp, 1, icomm, ISTATUS, ierr)
            zau(irct(2):irct(2) + irct(1) - 1) = zau(irct(2):irct(2) + irct(1) - 1) + wwr(:irct(1))
            wws(:irct(1)) = wwr(:irct(1))
            isct(:2) = irct(:2)
        enddo
    endif
    !$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp


!***HACApK_adot_lfmtx_hyp
subroutine HACApK_adot_lfmtx_hyp_32(zau, st_leafmtxp_32, st_ctl, zu_32, wws, wwr, isct, irct, nd)
    include 'mpif.h'
    type(st_HACApK_leafmtxp_32) :: st_leafmtxp_32
    type(st_HACApK_lcontrol) :: st_ctl
    real*8 :: zau(*), wws(*), wwr(*)
    real*4 :: zu_32(*)
    integer*4 :: isct(*), irct(*)
    integer*4 :: ISTATUS(MPI_STATUS_SIZE)
    integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
    real*8 :: time_start, time_end = 0.0
    ! integer*4, dimension(:), allocatable :: ISTATUS
    1000 format(5(a, i10)/)
    2000 format(5(a, f10.4)/)

    lpmd => st_ctl%lpmd(:)
    lnp(0:) => st_ctl%lnp
    lsp(0:) => st_ctl%lsp
    lthr(0:) => st_ctl%lthr
    ! allocate(ISTATUS(MPI_STATUS_SIZE))
    mpinr = lpmd(3)
    mpilog = lpmd(4)
    nrank = lpmd(2)
    icomm = lpmd(1)
    ndnr_s = lpmd(6)
    ndnr_e = lpmd(7)
    ndnr = lpmd(5)
    zau(:nd) = 0.0d0

    !$omp barrier

    !$omp master
      time_start = MPI_Wtime()
    !$omp end master

    do i = 1, 1000
      call HACApK_adot_body_lfmtx_hyp_rise_32(zau, st_leafmtxp_32, st_ctl, zu_32, nd)
    enddo

    !$omp barrier

    !$omp master
      time_end = MPI_Wtime()
      print *, "Total execution time:"
      print *, time_end - time_start
    !$omp end master

    !$omp barrier

    ! !$omp master
    ! print *, 'DEBUG 64'

    ! open(10, file='result_vector.data')
    !   do i = 1, nd
    !     write(10, *) zau(i)
    !   enddo
    ! close(10)

    ! !$omp end master
    ! !$omp barrier

    stop

    !$omp master
    if (nrank > 1) then
        wws(1:lnp(mpinr)) = zau(lsp(mpinr):lsp(mpinr) + lnp(mpinr) - 1)
        ncdp = mod(mpinr + 1, nrank)
        ncsp = mod(mpinr + nrank - 1, nrank)
        isct(1) = lnp(mpinr);isct(2) = lsp(mpinr)
        do ic = 1,nrank-1
            call MPI_SENDRECV(isct, 2, MPI_INTEGER, ncdp, 1, irct,2, MPI_INTEGER, ncsp, 1, icomm, ISTATUS, ierr)
            call MPI_SENDRECV(wws, isct, MPI_DOUBLE_PRECISION, ncdp, 1, wwr, irct, MPI_DOUBLE_PRECISION, ncsp, 1, icomm, ISTATUS, ierr)
            zau(irct(2):irct(2) + irct(1) - 1) = zau(irct(2):irct(2) + irct(1) - 1) + wwr(:irct(1))
            wws(:irct(1)) = wwr(:irct(1))
            isct(:2) = irct(:2)
        enddo
    endif
    !$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp_32


 !***HACApK_adot_lfmtx_hyp
subroutine HACApK_adot_lfmtx_hyp_32_64(zau, st_leafmtxp_32, st_ctl, zu, wws, wwr, isct, irct, nd)
    include 'mpif.h'
    type(st_HACApK_leafmtxp_32) :: st_leafmtxp_32
    type(st_HACApK_lcontrol) :: st_ctl
    real*8 :: zau(*), wws(*), wwr(*)
    real*8 :: zu(*)
    integer*4 :: isct(*), irct(*)
    integer*4 :: ISTATUS(MPI_STATUS_SIZE)
    integer*4, pointer :: lpmd(:), lnp(:), lsp(:), lthr(:)
    real*8 :: time_start, time_end = 0.0
    ! integer*4, dimension(:), allocatable :: ISTATUS
    1000 format(5(a, i10)/)
    2000 format(5(a, f10.4)/)

    lpmd => st_ctl%lpmd(:)
    lnp(0:) => st_ctl%lnp
    lsp(0:) => st_ctl%lsp
    lthr(0:) => st_ctl%lthr
    ! allocate(ISTATUS(MPI_STATUS_SIZE))
    mpinr = lpmd(3)
    mpilog = lpmd(4)
    nrank = lpmd(2)
    icomm = lpmd(1)
    ndnr_s = lpmd(6)
    ndnr_e = lpmd(7)
    ndnr = lpmd(5)
    zau(:nd) = 0.0d0

    !$omp barrier

    !$omp master
      time_start = MPI_Wtime()
    !$omp end master

    do i = 1, 1000
      call HACApK_adot_body_lfmtx_hyp_rise_32_64(zau, st_leafmtxp_32, st_ctl, zu, nd)
    enddo

    !$omp barrier

    !$omp master
      time_end = MPI_Wtime()
      print *, "Total execution time:"
      print *, time_end - time_start
    !$omp end master

    !$omp barrier

    ! !$omp master
    ! print *, 'DEBUG 64'

    ! open(10, file='result_vector.data')
    !   do i = 1, nd
    !     write(10, *) zau(i)
    !   enddo
    ! close(10)

    ! !$omp end master
    ! !$omp barrier

    stop

    !$omp master
    if (nrank > 1) then
        wws(1:lnp(mpinr)) = zau(lsp(mpinr):lsp(mpinr) + lnp(mpinr) - 1)
        ncdp = mod(mpinr + 1, nrank)
        ncsp = mod(mpinr + nrank - 1, nrank)
        isct(1) = lnp(mpinr);isct(2) = lsp(mpinr)
        do ic = 1,nrank-1
            call MPI_SENDRECV(isct, 2, MPI_INTEGER, ncdp, 1, irct,2, MPI_INTEGER, ncsp, 1, icomm, ISTATUS, ierr)
            call MPI_SENDRECV(wws, isct, MPI_DOUBLE_PRECISION, ncdp, 1, wwr, irct, MPI_DOUBLE_PRECISION, ncsp, 1, icomm, ISTATUS, ierr)
            zau(irct(2):irct(2) + irct(1) - 1) = zau(irct(2):irct(2) + irct(1) - 1) + wwr(:irct(1))
            wws(:irct(1)) = wwr(:irct(1))
            isct(:2) = irct(:2)
        enddo
    endif
    !$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp_32_64


!***HACApK_adot_lfmtx_hyp
    subroutine HACApK_adot_lfmtx_hyp_kawamura(zau, st_leafmtxp, st_ctl, zu, wws, wwr, isct, irct, nd)
        include 'mpif.h'
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl
        real*8 :: zau(*), zu(*), wws(*), wwr(*)
        integer*4 :: isct(*), irct(*)

        type(st_HACApK_leafmtxp), save :: newLeafmtxp_shared
        type(newLeafmtxpParam) :: NLMparam
        real*8 :: time_start, time_end, time_total = 0.0
        real*8 :: z1(110), z2(110), z3(110), z4(110), z5(110)
        integer*4 :: ISTATUS(MPI_STATUS_SIZE), values(8)
        integer*4, pointer :: lpmd(:), lnp(:), lsp(:)

        integer, parameter :: argc = 6, ALLMETHODS = 10, DIVISION = 15

        character(20), dimension(argc) :: argv
        character*20, dimension(ALLMETHODS) :: method_filename
        integer :: i, nThread, myThread, eof, fileA, fileB, readErr
        integer :: methodNo, chunksize, div
        integer :: divided, bigmatrix, rank1bigmatrix, rank2bigmatrix
        integer :: method_handler(ALLMETHODS)
        integer*4, allocatable :: threadStartIndex(:), threadEndIndex(:)
        real*8, allocatable, dimension(:) :: fileA_result, fileB_result

!   # args:: "mpiexec.hydra" - "./bem-bb-SCM.out" - –â‘èƒtƒ@ƒCƒ‹ - "output.lst" - ƒƒ\ƒbƒhNo - chunksize - div
!   # args:: "mpiexec.hydra" - 0th arg            - 1st arg      - 2nd arg      - 3rd arg    - 4th arg   - 5th arg
!   # args:: "mpiexec.hydra" - argv(1)            - argv(2)      - argv(3)      - argv(4)    - argv(5)   - argv(6)

        do i = 1, argc
            call getarg(i, argv(i + 1))
        enddo

        read(argv(4), *)    methodNo
        read(argv(5), *)    chunksize
        read(argv(6), *)    div

        nThread = omp_get_num_threads()
        myThread = omp_get_thread_num()
        eof = 0
        divided = 0;        bigmatrix = 0;      rank1bigmatrix = 0; rank2bigmatrix = 0

        zau(:nd) = 0.0d0
        zu(:nd) = 1d-6

!$omp single
        call makeNewLeafmtxp_test(st_leafmtxp, newLeafmtxp_shared, div)
        print *, "end makeNewLeafmtxp"
!$omp end single

        allocate(threadStartIndex(nThread), threadEndIndex(nThread))
        call methodFiveStaticOrderAllocation(newLeafmtxp_shared, nThread, threadStartIndex, threadEndIndex) ! method9

        if(methodNo == 2 .or. methodNo == 4) then
            do ip = 0, nd
                if(st_leafmtxp%st_lf(ip)%kt .gt. DIVISION) then
                    bigmatrix = bigmatrix + 1
                    if(st_leafmtxp%st_lf(ip)%ltmtx == 1)    rank1bigmatrix = rank1bigmatrix + 1
                    if(st_leafmtxp%st_lf(ip)%ltmtx == 2)    rank2bigmatrix = rank2bigmatrix + 1
                endif
                if (st_leafmtxp%st_lf(ip)%ltmtx==1 .and. (st_leafmtxp%st_lf(ip)%kt .gt. DIVISION)) then
                    st_leafmtxp%st_lf(ip)%ltmtx = 3
                    divided = divided + 1
                endif
            enddo
        endif

!$omp single
        print *, "end allocation"
        print *, ""
        print *, "problem file   ", argv(2)
        if(methodNo == 2 .or. methodNo ==4) then
            print *, "method 2 4 division  ", DIVISION
            print *, "method 2 4 dividedBigmatrix ", divided, " bigmatrix ", bigmatrix
            print *, "bigmatrix rank 1, 2       ", rank1bigmatrix, rank2bigmatrix
        endif
        print '("threads ", i2, "  methodNo ", i2, "  chunksize ", i5, "  div ", i5)',  nThread, methodNo, chunksize, div
        print *, ""
!$omp end single


        time_start = MPI_Wtime()
!$omp barrier
        ! do i = 1, 1
        do i = 1, 1
            select case(methodNo)
                case(1)
                    call HACApK_adot_body_lfmtx_hyp_original(zau, st_leafmtxp, st_ctl, zu, nd)              ! method1
                case(2)
                    call HACApK_adot_body_lfmtx_hyp(zau, st_leafmtxp, st_ctl, zu, nd)                       ! method2
                case(3)
                    call HACApK_adot_body_lfmtx_hyp_dynamic(zau, st_leafmtxp, st_ctl, zu, nd, chunksize)    ! method3
                case(4)
                    call HACApK_adot_body_lfmtx_hyp2(zau, st_leafmtxp, st_ctl, zu, nd, chunksize)           ! method4
                case(5)
                    call HACApK_adot_body_lfmtx_expand(zau, newLeafmtxp_shared, zu, nd, chunksize)          ! method5
                case(6)
                    call HACApK_adot_body_lfmtx_expandB2(zau, newLeafmtxp_shared, zu, nd, chunksize, threadStartIndex(myThread + 1), threadEndIndex(myThread + 1))  ! method9 = method5B-2
                case(9)
                    call HACApK_adot_body_lfmtx_expandB2(zau, newLeafmtxp_shared, zu, nd, chunksize, threadStartIndex(myThread + 1), threadEndIndex(myThread + 1))  ! method9 = method5B-2
                ! case(24)
                !     call HACApK_adot_body_lfmtx_hyp_rise_32(zau, st_leafmtxp_32, st_ctl, zu_32, nd)
                case(25)
                    call HACApK_adot_body_lfmtx_hyp_rise_64(zau, st_leafmtxp, st_ctl, zu, nd)
                case default
                    if(myThread==0)     print *, "no methods done"
            end select
        enddo
        time_end = MPI_Wtime()
        time_total = time_end - time_start
!$omp barrier

        deallocate(threadStartIndex, threadEndIndex)

!$OMP MASTER

                ! open(10, file='result_vector.data')

                !      do i = 1, nd
                !        write(10, *) zau(i)
                !      enddo

                ! close(10)

        print '("total execution time:    ",  f11.6)', time_total
        print *, ""

        do i = 1, ALLMETHODS
            write (method_filename(i), '("method", i0, "_", i0, ".dat")')    i, nd
            method_handler(i) = 10 + i
        enddo

        ! method‚ÌŒ‹‰Ê‚ð‘‚«ž‚Ý
        if(.false.) then
            open(method_handler(methodNo), file = method_filename(methodNo), status = "replace")
            do i = 1, st_leafmtxp%nlf
                write(method_handler(methodNo), *), zau(i)
            enddo
            close(method_handler(methodNo))
            print *, "wrote ", method_filename(methodNo)
            print *, ""
        endif

        ! method‚²‚Æ‚ÌŒvŽZŒ‹‰Êƒtƒ@ƒCƒ‹”äŠr
        if(.false.) then
            fileA = 1
            fileB = methodNo

            print '("START ", a20, " and ", a20, " difference check")', method_filename(fileA), method_filename(fileB)

            open(method_handler(fileA), file = method_filename(fileA), status = "old", IOSTAT = readErr)
            open(method_handler(fileB), file = method_filename(fileB), status = "old", IOSTAT = readErr)

            if(readErr == 0) then
        ! method‚ÌŒ‹‰Ê‚ð“Ç‚Ýž‚Ý
                allocate(fileA_result(st_leafmtxp%nlf), fileB_result(st_leafmtxp%nlf))
                do i = 1, st_leafmtxp%nlf
                    read(method_handler(fileA), *, IOSTAT = eof), fileA_result(i)
                    read(method_handler(fileB), *, IOSTAT = eof), fileB_result(i)
                    if(eof .ne. 0)  exit
                enddo

        ! method‚ÌŒ‹‰Ê‚ð”äŠr‚·‚éi‚ ‚Ü‚è‚É‚à”’l‚É·‚ª‚ ‚éê‡‚ÍŒvŽZŒ‹‰Ê‚ªŠÔˆá‚Á‚Ä‚¢‚éj
                do i = 1, st_leafmtxp%nlf
                    if(fileB_result(i)==0) then                 !!!! do nothing
                    elseif( abs(fileA_result(i) - fileB_result(i)) / abs(fileB_result(i)) > 10d-6) then             !!!! Œë·è‡’l
                        print *, "differs at index ", i, " ", abs(fileA_result(i) - fileB_result(i)) / abs(fileB_result(i))
                    endif
                enddo
                print '("END ", a20, " and ", a20, " difference check", /)', method_filename(fileA), method_filename(fileB)
                deallocate(fileA_result, fileB_result)
            else
                print '("file reading error", /)'
            endif
            close(method_handler(fileA))
            close(method_handler(fileB))
        endif

        stop 'file difference check'
!$OMP END MASTER

!$omp master
        lpmd => st_ctl%lpmd(:)
        lnp(0:) => st_ctl%lnp
        lsp(0:) => st_ctl%lsp
        mpinr = lpmd(3)
        nrank = lpmd(2)
        icomm = lpmd(1)

        if(nrank>1) then
            wws(1:lnp(mpinr)) = zau(lsp(mpinr):lsp(mpinr) + lnp(mpinr)-1)
            ncdp = mod(mpinr + 1, nrank)
            ncsp = mod(mpinr + nrank-1, nrank)
            isct(1) = lnp(mpinr)
            isct(2) = lsp(mpinr)

            do ic = 1, nrank-1
                call MPI_SENDRECV(isct, 2,MPI_INTEGER, ncdp, 1, irct, 2, MPI_INTEGER, ncsp, 1, icomm, ISTATUS, ierr)
                call MPI_SENDRECV(wws, isct, MPI_DOUBLE_PRECISION, ncdp, 1, wwr, irct, MPI_DOUBLE_PRECISION, ncsp, 1, icomm, ISTATUS, ierr)
                zau(irct(2):irct(2) + irct(1) - 1) = zau(irct(2):irct(2) + irct(1)-1) + wwr(:irct(1))
                wws(:irct(1)) = wwr(:irct(1))
                isct(:2) = irct(:2)
            enddo
        endif
!$omp end master
    end subroutine HACApK_adot_lfmtx_hyp_kawamura


!***HACApK_adot_body_lfmtx_hyp_setLtmtx
    subroutine HACApK_adot_body_lfmtx_hyp_setLtmtx(st_leafmtxp, st_ctl)
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl

        nThread = omp_get_num_threads()
        ipFirst = st_ctl%lthr(0)
        ipLast  = st_ctl%lthr(nThread) - 1

        do ip = ipFirst, ipLast
            if (st_leafmtxp%st_lf(ip)%ltmtx==1 .and. (st_leafmtxp%st_lf(ip)%kt .gt. nThread*6)) then
                st_leafmtxp%st_lf(ip)%ltmtx = 3
            endif
        enddo
    endsubroutine HACApK_adot_body_lfmtx_hyp_setLtmtx


! include "m_HACApK_solve_methods.f90"

!***HACApK_adot_pmt_lfmtx_hyp
    integer function HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp, st_bemv, st_ctl, aww, ww)
        include 'mpif.h'
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl
        type(st_HACApK_calc_entry) :: st_bemv
        real*8 :: ww(*), aww(*)
        real*8, dimension(:), allocatable :: u, au, wws, wwr
        integer*4, dimension(:), allocatable :: isct, irct
        integer*4, pointer :: lnp(:)

        lrtrn = 0
        lnp(0:) => st_ctl%lnp
        nd = st_bemv%nd
        icomm = st_ctl%lpmd(1)
        nrank = st_ctl%lpmd(2)

        allocate(u(nd), au(nd), isct(2), irct(2))
        u(:nd) = ww(st_ctl%lod(:nd))
        allocate(wws(maxval(st_ctl%lnp(:nrank))), wwr(maxval(st_ctl%lnp(:nrank))))
        call MPI_Barrier(icomm, ierr)
!$omp parallel
!$omp barrier
        call HACApK_adot_lfmtx_hyp(au, st_leafmtxp, st_ctl, u, wws, wwr, isct, irct, nd)
!$omp barrier
!$omp end parallel
        call MPI_Barrier( icomm, ierr )
        aww(st_ctl%lod(:nd)) = au(:nd)
        HACApK_adot_pmt_lfmtx_hyp = lrtrn
    end function HACApK_adot_pmt_lfmtx_hyp


!***HACApK_adot_pmt_lfmtx_p
    integer function HACApK_adot_pmt_lfmtx_p(st_leafmtxp, st_bemv, st_ctl, aww, ww)
        include 'mpif.h'
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl
        type(st_HACApK_calc_entry) :: st_bemv
        real*8 :: ww(st_bemv%nd), aww(st_bemv%nd)
        real*8, dimension(:), allocatable :: u, au

        lrtrn = 0
        icomm = st_ctl%lpmd(1)
        nd = st_bemv%nd

        allocate(u(nd), au(nd))
        u(:nd) = ww(st_ctl%lod(:nd))

        call MPI_Barrier(icomm, ierr)
        call HACApK_adot_lfmtx_p(au, st_leafmtxp, st_ctl, u, nd)
        aww(st_ctl%lod(:nd)) = au(:nd)
        HACApK_adot_pmt_lfmtx_p = lrtrn
    end function HACApK_adot_pmt_lfmtx_p


!***HACApK_measurez_time_ax_lfmtx
    subroutine HACApK_measurez_time_ax_lfmtx(st_leafmtxp, st_ctl, nd, nstp, lrtrn)
        include 'mpif.h'
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_lcontrol) :: st_ctl

        real*8, dimension(:), allocatable :: wws, wwr, u, b
        integer*4 :: isct(2), irct(2)
        integer*4, pointer :: lnp(:)

        lnp(0:) => st_ctl%lnp
        nrank = st_ctl%lpmd(2)
        mstep = st_ctl%param(99)

        allocate(u(nd), b(nd), wws(maxval(lnp(0:nrank-1))), wwr(maxval(lnp(0:nrank-1))))
!$omp parallel private(il)
        do il = 1, mstep
            u(:) = 1.0
            b(:) = 1.0
            call HACApK_adot_lfmtx_hyp(u, st_leafmtxp, st_ctl, b, wws, wwr, isct, irct, nd)
        enddo
!$omp end parallel
        deallocate(wws, wwr)
    end subroutine HACApK_measurez_time_ax_lfmtx

endmodule m_HACApK_solve

