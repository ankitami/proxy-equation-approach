!*****************************************************************************************************************************************************!
! 3D Advection Diffusion Equation : Asynchronous Implementation (Proxy-Equation approach) 
! MPI Code to compare the following
!
! 1. Async Governing Equation Scheme
! 2. Synchronous Scheme
!
! DATE LAST MODIFIED : July 27, 2017
!
! AUTHOR : Ankita Mittal
!*****************************************************************************************************************************************************!
program main
        implicit none
        include "mpif.h"
        real*8, parameter :: norm_t = 1.0D0
        real*8, parameter :: pi = 3.14159265358979323846264338D0
        real*8, parameter :: alpha = 0.1D0
        real*8, parameter :: cx = 1.0D0
        real*8, parameter :: cy = 0.0D0
        real*8, parameter :: cz = 0.0D0
        real*8, parameter :: r_alpha = 0.1D0

        real*8 :: angle
        integer*4 :: K_x, K_y, K_z
        integer*4, parameter :: K_num = 1
        integer*4, parameter :: angle_num = 1

        integer*4, parameter :: L = 2
        integer*4, parameter :: flag = 0

        integer*4, parameter :: prob_num = 2
                                                !     X,     Y,     Z    Directional Delay Probabilities
        real*8, dimension (3,prob_num) :: prob0 = (/1.0D0, 1.0D0, 1.0D0, &
                                                    0.6D0, 0.6D0, 0.6D0/)!, &
                                                  !  0.3D0, 0.3D0, 0.3D0, &
                                                  !  0.0D0, 0.0D0, 0.0D0/)

        integer*4, parameter :: n_num = 3
        integer*4, dimension (n_num) :: nx = (/32, 64, 128/)!, 256/)!, 512/)!, 1024/)
        integer*4, dimension (n_num) :: ny = (/32, 64, 128/)!, 256/)!, 512/)!, 1024/)
        integer*4, dimension (n_num) :: nz = (/32, 64, 128/)!, 256/)!, 512/)!, 1024/)

        integer*4, parameter :: n_write_sol = 4
        real*8, dimension(n_write_sol) :: write_normt = (/0.25D0,0.5D0,0.75D0,1.0D0/)

        real*8, dimension (3,0:L-1) :: prob, prob_cum
        real*8, dimension(K_num,n_num) :: Avg_aerr, Avg_err, Avg_all_aerr, Avg_all_err
        real*8, dimension(n_num) :: ref1, ref2
        real*8 :: r_cx, r_cy, r_cz, alpha_new, cx_new, cy_new, cz_new, k_fac
	real*8 :: dx, dy, dz, dt

	logical :: periods(3), reorder
        real*8, allocatable :: u_exact(:,:,:), x(:), y(:), z(:)
        real*8, allocatable :: u_aproc(:,:,:,:), u_proc(:,:,:,:)
        real*8, allocatable :: u_aerr(:,:,:,:), u_err(:,:,:,:)
        real*8, allocatable :: async_d(:,:,:)
        real*8 :: rand(3), time
        real*8, dimension(angle_num) :: rand_angle, phis

        CHARACTER (len=50) :: filename
        CHARACTER(*), PARAMETER :: fileplace = "./"

        integer*4 :: n_per_px, n_per_py, n_per_pz
        integer*4 :: error, nproc, taskid, rank_cart, COMM_CART
        integer*4 :: status(12), dims(3), rextent, rextent1
        integer*4 :: i, j, k, n, ln, steps_proc, prob_par, nprocs, T, angle_par, k_par, jj, ii, i_write, t_par
        integer*4 :: nprocx, nprocy, nprocz, coords(3), cord_neh(3)
        integer*4 :: left_right, back_front, down_up, temp1
        integer*4 :: left, right, up, down, front, back
        integer*4 :: k_l, k_r, k_u, k_d, k_f, k_b, k_del(3)
        integer*4 :: req_ap(12), req_p(12), req_k(12)

        ! Initialize MPI
        call MPI_Init(error)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank_cart, error)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, error)

        dims(1) = 0
        dims(2) = 0
        dims(3) = 0

        call MPI_DIMS_CREATE(nproc, 3, dims, error)

        if(rank_cart == 0) then
            write(*,*) 'rank_cart',rank_cart, 'nproc', nproc, 'Dims', dims(1), dims(2), dims(3)
        endif
        periods(1) = .TRUE.
        periods(2) = .TRUE.
        periods(3) = .TRUE.
        reorder = .TRUE.

        call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, periods, reorder, COMM_CART, error)
        call MPI_COMM_RANK(COMM_CART, taskid, error)
        call MPI_CART_COORDS(COMM_CART, taskid, 3, coords, error)

        do angle_par = 1,angle_num
            phis(angle_par) = -16.5D0 !2.0D1*(rand_angle(angle_par)-0.5D0)*2.0D0
        enddo

        do i = 1,12
            req_ap(i) = i
            req_p(i) = i+12
            req_k(i) = i+24
        enddo

        do prob_par = 1, prob_num
            do k = 1,3
                prob(k,0) = prob0(k,prob_par)
                prob(k,1) = 1.0D0 - prob0(k,prob_par)
                prob_cum(k,0) = prob(k,0)
                do n = 1, L-1
                    prob_cum(k,n) = prob_cum(k,n-1) + prob(k,n)
                enddo
            enddo
                if(taskid == 0) then
                        write ( *, * ) ' Probabilities at each delay = ', prob
                        write ( *, * ) ' Cumulative Probabilities = ', prob_cum
                endif
        do n = 1,n_num
                dx = 2.0D0*pi/nx(n)
                dy = 2.0D0*pi/ny(n)
                dz = 2.0D0*pi/nz(n)
                k_fac = 2.0D0
                nprocx = dims(1)
                nprocy = dims(2)
                nprocz = dims(3)
                n_per_px = nx(n)/nprocx
                n_per_py = ny(n)/nprocy
                n_per_pz = nz(n)/nprocz

                dt = min(r_alpha*dx*dx/alpha, r_alpha*dy*dy/alpha, r_alpha*dz*dz/alpha)
                r_cx = cx*dt/dx
                r_cy = cy*dt/dy
                r_cz = cz*dt/dz

                T = INT(FLOOR(norm_t*2.0D0*pi/(dabs(max(cx,cy,cz))*dt)))
                if(taskid == 0) then
                    write(*,*) 'T = ', T, 'nx = ', nx(n), 'ny = ', ny(n), 'nz = ', nz(n)
                endif

                allocate(u_aproc(0:n_per_px+1,0:n_per_py+1,0:n_per_pz+1,0:L))
                allocate(u_aerr(0:n_per_px+1,0:n_per_py+1,0:n_per_pz+1,0:L))
                allocate(u_proc(0:n_per_px+1,0:n_per_py+1,0:n_per_pz+1,0:L))
                allocate(u_err(0:n_per_px+1,0:n_per_py+1,0:n_per_pz+1,0:L))
                allocate(async_d(n_per_px,n_per_py,n_per_pz))

                allocate(u_exact(n_per_px,n_per_py,n_per_pz))
                allocate(x(n_per_px))
                allocate(y(n_per_py))
                allocate(z(n_per_pz))

                !************************************************************************!
                ! Communication Data types:
                ! 1. back_front ==> xy plane
                ! 2. down_up    ==> xz plane
                ! 3. left_right ==> yz plane
                !************************************************************************!

                ! Define type for front-back transfer
                call MPI_TYPE_VECTOR(n_per_py, n_per_px, n_per_px+2, MPI_REAL8, back_front, error)
                call MPI_TYPE_COMMIT(back_front, error)

                ! Define type for up-down transfer
                call MPI_TYPE_VECTOR(n_per_pz, n_per_px, (n_per_px+2)*(n_per_py+2), MPI_REAL8, down_up, error)
                call MPI_TYPE_COMMIT(down_up, error)

                ! Define type for left-right transfer
                call MPI_TYPE_VECTOR(n_per_py, 1, n_per_px+2, MPI_REAL8, temp1, error)
                call MPI_TYPE_COMMIT(temp1, error)
                call MPI_TYPE_EXTENT(MPI_REAL8, rextent, error)
                call MPI_TYPE_EXTENT(temp1, rextent1, error)
                call MPI_TYPE_HVECTOR(n_per_pz, 1, (n_per_px+2)*(n_per_py+2)*rextent, temp1, left_right, error)
                call MPI_TYPE_COMMIT(left_right, error)

                write(*,*) 'taskid', taskid, 'Coords', coords(1), coords(2), coords(3)

        do k_par = 1,K_num
        
                K_x = k_par*k_fac
                K_y = k_par*k_fac
                K_z = k_par*k_fac

                if(taskid == 0) then
                        write ( *, * ) ' Wavenumbers (kx,ky,kz) = ', K_x, K_y, K_z
                endif

                Avg_aerr(k_par,n) = 0.0D0
                Avg_err(k_par,n) = 0.0D0
                Avg_all_aerr(k_par,n) = 0.0D0
                Avg_all_err(k_par,n) = 0.0D0
                ref1(n) = 0.0D0
                ref2(n) = 0.0D0

        do angle_par = 1,angle_num

                angle = phis(angle_par)

!*****************************************************************************************************************************************************!
                ! Initialization
!*****************************************************************************************************************************************************!
                do k = 1, n_per_pz
                        z(k) = (coords(3)*n_per_pz+k-1)*dz
                        do j = 1, n_per_py
                                y(j) = (coords(2)*n_per_py+j-1)*dy
                                do i = 1, n_per_px
                                        x(i) = (coords(1)*n_per_px+i-1)*dx
                                        u_exact(i,j,k) = dsin(K_x*x(i)+K_y*y(j)+K_z*z(k)+pi*angle/180.0D0)

                                        u_aproc(i,j,k,0) = u_exact(i,j,k)
                                        u_aerr(i,j,k,0) = u_aproc(i,j,k,0) - u_exact(i,j,k)

                                        u_proc(i,j,k,0) = u_exact(i,j,k)
                                        u_err(i,j,k,0) = u_proc(i,j,k,0) - u_exact(i,j,k)

                                        async_d(i,j,k) = 0.0D0
                                enddo
                        enddo
                enddo

                !write(*,*) 'taskid', taskid, 'x = ', x
                !write(*,*) 'taskid', taskid, 'y = ', y
                !write(*,*) 'taskid', taskid, 'z = ', z

                cord_neh(1) = coords(1)
                cord_neh(2) = coords(2)
                cord_neh(3) = mod(coords(3)+1,nprocz)
                call MPI_CART_RANK(COMM_CART, cord_neh, front, error)

                cord_neh(1) = coords(1)
                cord_neh(2) = coords(2)
                cord_neh(3) = mod(coords(3)-1,nprocz)
                call MPI_CART_RANK(COMM_CART, cord_neh, back, error)

                cord_neh(1) = coords(1)
                cord_neh(2) = mod(coords(2)+1,nprocy)
                cord_neh(3) = coords(3)
                call MPI_CART_RANK(COMM_CART, cord_neh, up, error)

                cord_neh(1) = coords(1)
                cord_neh(2) = mod(coords(2)-1,nprocy)
                cord_neh(3) = coords(3)
                call MPI_CART_RANK(COMM_CART, cord_neh, down, error)

                cord_neh(1) = mod(coords(1)+1,nprocx)
                cord_neh(2) = coords(2)
                cord_neh(3) = coords(3)
                call MPI_CART_RANK(COMM_CART, cord_neh, right, error)

                cord_neh(1) = mod(coords(1)-1,nprocx)
                cord_neh(2) = coords(2)
                cord_neh(3) = coords(3)
                call MPI_CART_RANK(COMM_CART, cord_neh, left, error)

                !if(taskid == 0) then
                !    write(*,*) 'taskid', taskid, 'right = ', right, 'Left = ', left
                !    write(*,*) 'taskid', taskid, 'up = ', up, 'down = ', down
                !    write(*,*) 'taskid', taskid, 'front = ', front, 'back = ', back
                !endif

                time = 0.0D0

                call MPI_BARRIER(MPI_COMM_WORLD, error)
                
        ! Loop over all time-steps
        do t_par = 0,T-1
                time = time + dt
                call Random_number(rand)
                
!*****************************************************************************************************************************************************!
                ! Delay Levels Set up
!*****************************************************************************************************************************************************!
                do k = 1,3
                    k_del(k) = 0
                    do ln = 1,L-1
                        if((rand(k) > prob_cum(k,ln-1)) .AND. (rand(k) <= prob_cum(k,ln))) THEN
                            k_del(k) = min(t_par,ln)
                        endif
                    enddo
                enddo
                !write(*,*) 'taskid = ', taskid, 'k_del = ', k_del

!*****************************************************************************************************************************************************!
                ! Exact Solution
!*****************************************************************************************************************************************************!
                do k = 1,n_per_pz
                    do j = 1,n_per_py
                        do i = 1,n_per_px
                            u_exact(i,j,k) = DEXP(-time*alpha*(K_x*K_x+K_y*K_y+K_z*K_z))*DSIN(K_x*(x(i)-cx*time) + &
                                                     K_y*(y(k)-cy*time) + K_z*(z(k)-cz*time) + pi*angle/180.0D0)
                        enddo
                    enddo
                enddo

!*****************************************************************************************************************************************************!
                ! Communications
!*****************************************************************************************************************************************************!
                ! Corrected Data
                call MPI_IRECV(u_aproc(1,1,0,mod(t_par,L+1)), 1, back_front, back, 0, COMM_CART, req_ap(1), error)
                call MPI_ISEND(u_aproc(1,1,1,mod(t_par-k_del(3),L+1)), 1, back_front, back, 1, COMM_CART, req_ap(2), error)
                call MPI_IRECV(u_aproc(1,1,n_per_pz+1,mod(t_par,L+1)), 1, back_front, front, 1, COMM_CART, req_ap(3), error)
                call MPI_ISEND(u_aproc(1,1,n_per_pz,mod(t_par-k_del(3),L+1)), 1, back_front, front, 0, COMM_CART, req_ap(4), error)
                call MPI_IRECV(u_aproc(1,0,1,mod(t_par,L+1)), 1, down_up, down, 3, COMM_CART, req_ap(5), error)
                call MPI_ISEND(u_aproc(1,1,1,mod(t_par-k_del(2),L+1)), 1, down_up, down, 2, COMM_CART, req_ap(6), error)
                call MPI_IRECV(u_aproc(1,n_per_py+1,1,mod(t_par,L+1)), 1, down_up, up, 2, COMM_CART, req_ap(7), error)
                call MPI_ISEND(u_aproc(1,n_per_py,1,mod(t_par-k_del(2),L+1)), 1, down_up, up, 3, COMM_CART, req_ap(8), error)
                call MPI_IRECV(u_aproc(0,1,1,mod(t_par,L+1)), 1, left_right, left, 5, COMM_CART, req_ap(9), error)
                call MPI_ISEND(u_aproc(1,1,1,mod(t_par-k_del(1),L+1)), 1, left_right, left, 4, COMM_CART, req_ap(10), error)
                call MPI_IRECV(u_aproc(n_per_px+1,1,1,mod(t_par,L+1)), 1, left_right, right, 4, COMM_CART, req_ap(11), error)                
                call MPI_ISEND(u_aproc(n_per_px,1,1,mod(t_par-k_del(1),L+1)), 1, left_right, right, 5, COMM_CART, req_ap(12), error)

                ! Delay Data
                call MPI_IRECV(k_b, 1, MPI_INTEGER4, back, 7, COMM_CART, req_k(1), error)
                call MPI_ISEND(k_del(3), 1, MPI_INTEGER4, back, 8, COMM_CART, req_k(2), error)
                call MPI_IRECV(k_f, 1, MPI_INTEGER4, front, 8, COMM_CART, req_k(3), error)
                call MPI_ISEND(k_del(3), 1, MPI_INTEGER4, front, 7, COMM_CART, req_k(4), error)
                call MPI_IRECV(k_d, 1, MPI_INTEGER4, down, 9, COMM_CART, req_k(5), error)
                call MPI_ISEND(k_del(2), 1, MPI_INTEGER4, down, 10, COMM_CART, req_k(6), error)
                call MPI_IRECV(k_u, 1, MPI_INTEGER4, up, 10, COMM_CART, req_k(7), error)
                call MPI_ISEND(k_del(2), 1, MPI_INTEGER4, up, 9, COMM_CART, req_k(8), error)
                call MPI_IRECV(k_l, 1, MPI_INTEGER4, left, 11, COMM_CART, req_k(9), error)
                call MPI_ISEND(k_del(1), 1, MPI_INTEGER4, left, 12, COMM_CART, req_k(10), error)
                call MPI_IRECV(k_r, 1, MPI_INTEGER4, right, 12, COMM_CART, req_k(11), error)
                Call MPI_ISEND(k_del(1), 1, MPI_INTEGER4, right, 11, COMM_CART, req_k(12), error)

                ! Non-Corrected Data
                call MPI_ISEND(u_proc(1,1,1,mod(t_par-k_del(3),L+1)), 1, back_front, back, 13, COMM_CART, req_p(1), error)
                call MPI_IRECV(u_proc(1,1,0,mod(t_par,L+1)), 1, back_front, back, 14, COMM_CART, req_p(2), error)
                call MPI_ISEND(u_proc(1,1,n_per_pz,mod(t_par-k_del(3),L+1)), 1, back_front, front, 14, COMM_CART, req_p(3), error)
                call MPI_IRECV(u_proc(1,1,n_per_pz+1,mod(t_par,L+1)), 1, back_front, front, 13, COMM_CART, req_p(4), error)
                call MPI_ISEND(u_proc(1,1,1,mod(t_par-k_del(2),L+1)), 1, down_up, down, 15, COMM_CART, req_p(5), error)
                call MPI_IRECV(u_proc(1,0,1,mod(t_par,L+1)), 1, down_up, down, 16, COMM_CART, req_p(6), error)
                call MPI_ISEND(u_proc(1,n_per_py,1,mod(t_par-k_del(2),L+1)), 1, down_up, up, 16, COMM_CART, req_p(7), error)
                call MPI_IRECV(u_proc(1,n_per_py+1,1,mod(t_par,L+1)), 1, down_up, up, 15, COMM_CART, req_p(8), error)
                call MPI_ISEND(u_proc(1,1,1,mod(t_par-k_del(1),L+1)), 1, left_right, left, 17, COMM_CART, req_p(9), error)
                call MPI_IRECV(u_proc(0,1,1,mod(t_par,L+1)), 1, left_right, left, 18, COMM_CART, req_p(10), error)
                call MPI_ISEND(u_proc(n_per_px,1,1,mod(t_par-k_del(1),L+1)), 1, left_right, right, 18, COMM_CART, req_p(11), error)
                call MPI_IRECV(u_proc(n_per_px+1,1,1,mod(t_par,L+1)), 1, left_right, right, 17, COMM_CART, req_p(12), error)

                call MPI_WaitALL(12, req_ap, status, error)
                call MPI_WaitALL(12, req_p, status, error)
                call MPI_WaitALL(12, req_k, status, error)

                !write(*,*) 'taskid = ', taskid, 'k recv = ', k_r, k_l, k_u, k_d, k_f, k_b

                !if(taskid == 0) then
                !        !write(*,*) 'taskid', taskid, 'all u', u_aproc
                !        write(*,*) 'taskid', taskid, 'x-delay', k_del(3)!u_aproc(1:n_per_px,1:n_per_py,n_per_pz,mod(t_par-k_del(1),L+1))
                !        write(*,*) 'taskid', taskid, 'y-delay', k_del(2)!u_aproc(1:n_per_px,1:n_per_py,1,mod(t_par-k_del(1),L+1))
                !        write(*,*) 'taskid', taskid, 'z-delay', k_del(1)!u_aproc(1:n_per_px,1,1:n_per_pz,mod(t_par-k_del(2),L+1))
                !        !write(*,*) 'taskid', taskid, 'up sent data', u_aproc(1:n_per_px,n_per_py,1:n_per_pz,mod(t_par-k_del(2),L+1))
                !        !write(*,*) 'taskid', taskid, 'left sent data', u_aproc(1,1:n_per_py,1:n_per_pz,mod(t_par-k_del(3),L+1))
                !        !write(*,*) 'taskid', taskid, 'right sent data', u_aproc(n_per_px,1:n_per_py,1:n_per_pz,mod(t_par-k_del(3),L+1))
                !endif

                !call MPI_BARRIER(MPI_COMM_WORLD, error)

                !if(taskid == 1) then
                !        !write(*,*) 'taskid', taskid, 'all u', u_aproc
                !        write(*,*) 'taskid', taskid, 'front delay', k_f!u_aproc(1:n_per_px,1:n_per_py,0,mod(t_par,L+1))
                !endif
                !if(taskid == 7) then
                !        write(*,*) 'taskid', taskid, 'back delay', k_b!u_aproc(1:n_per_px,1:n_per_py,n_per_pz+1,mod(t_par,L+1))
                !endif

                !call MPI_BARRIER(MPI_COMM_WORLD, error)

                !if(taskid == 8) then
                !        !write(*,*) 'taskid', taskid, 'all u', u_aproc
                !        write(*,*) 'taskid', taskid, 'up delay', k_u!u_aproc(1:n_per_px,0,1:n_per_pz,mod(t_par,L+1))
                !endif
                !if(taskid == 56) then
                !        write(*,*) 'taskid', taskid, 'down delay', k_d!u_aproc(1:n_per_px,n_per_py+1,1:n_per_pz,mod(t_par,L+1))
                !endif

                !call MPI_BARRIER(MPI_COMM_WORLD, error)

                !if(taskid == 64) then
                !        !write(*,*) 'taskid', taskid, 'all u', u_aproc
                !        write(*,*) 'taskid', taskid, 'right delay', k_r!u_aproc(0,1:n_per_py,1:n_per_pz,mod(t_par,L+1))
                !endif
                !if(taskid == 448) then
                !        write(*,*) 'taskid', taskid, 'left delay', k_l!u_aproc(n_per_px+1,1:n_per_py,1:n_per_pz,mod(t_par,L+1))
                !endif
                !call MPI_BARRIER(MPI_COMM_WORLD, error)

!*****************************************************************************************************************************************************!
                ! Commputations
!*****************************************************************************************************************************************************!                

                !************************************************************************!
                ! Determining Asynchronous Delay factor (async_d) throughout the domain
                !************************************************************************!

                ! Corners (8 points)
                async_d(1,1,1) = k_l*(r_alpha+r_cx/2.0D0)+k_d*(r_alpha+r_cy/2.0D0)+k_b*(r_alpha+r_cz/2.0D0)
                async_d(n_per_px,1,1) = k_r*(r_alpha-r_cx/2.0D0)+k_d*(r_alpha+r_cy/2.0D0)+k_b*(r_alpha+r_cz/2.0D0)
                async_d(1,n_per_py,1) = k_l*(r_alpha+r_cx/2.0D0)+k_u*(r_alpha-r_cy/2.0D0)+k_b*(r_alpha+r_cz/2.0D0)
                async_d(1,1,n_per_pz) = k_l*(r_alpha+r_cx/2.0D0)+k_d*(r_alpha+r_cy/2.0D0)+k_f*(r_alpha-r_cz/2.0D0)
                async_d(n_per_px,n_per_py,1) = k_r*(r_alpha-r_cx/2.0D0)+k_u*(r_alpha-r_cy/2.0D0)+k_b*(r_alpha+r_cz/2.0D0)
                async_d(n_per_px,1,n_per_pz) = k_r*(r_alpha-r_cx/2.0D0)+k_d*(r_alpha+r_cy/2.0D0)+k_f*(r_alpha-r_cz/2.0D0)
                async_d(1,n_per_py,n_per_pz) = k_l*(r_alpha+r_cx/2.0D0)+k_u*(r_alpha-r_cy/2.0D0)+k_f*(r_alpha-r_cz/2.0D0)
                async_d(n_per_px,n_per_py,n_per_pz) = k_r*(r_alpha-r_cx/2.0D0)+k_u*(r_alpha-r_cy/2.0D0)+k_f*(r_alpha-r_cz/2.0D0)

                ! Edges (12 lines)
                do i = 2, n_per_px-1
                        async_d(i,1,1) = k_d*(r_alpha+r_cy/2.0D0)+k_b*(r_alpha+r_cz/2.0D0)
                        async_d(i,n_per_py,1) = k_u*(r_alpha-r_cy/2.0D0)+k_b*(r_alpha+r_cz/2.0D0)
                        async_d(i,1,n_per_pz) = k_d*(r_alpha+r_cy/2.0D0)+k_f*(r_alpha-r_cz/2.0D0)
                        async_d(i,n_per_py,n_per_pz) = k_u*(r_alpha-r_cy/2.0D0)+k_f*(r_alpha-r_cz/2.0D0)
                enddo
                do j = 2, n_per_py-1
                        async_d(1,j,1) = k_l*(r_alpha+r_cx/2.0D0)+k_b*(r_alpha+r_cz/2.0D0)
                        async_d(n_per_px,j,1) = k_r*(r_alpha-r_cx/2.0D0)+k_b*(r_alpha+r_cz/2.0D0)
                        async_d(1,j,n_per_pz) = k_l*(r_alpha+r_cx/2.0D0)+k_f*(r_alpha-r_cz/2.0D0)
                        async_d(n_per_px,j,n_per_pz) = k_r*(r_alpha-r_cx/2.0D0)+k_f*(r_alpha-r_cz/2.0D0)
                enddo
                do k = 2, n_per_pz-1
                        async_d(1,1,k) = k_l*(r_alpha+r_cx/2.0D0)+k_d*(r_alpha+r_cy/2.0D0)
                        async_d(n_per_px,1,k) = k_r*(r_alpha-r_cx/2.0D0)+k_d*(r_alpha+r_cy/2.0D0)
                        async_d(1,n_per_py,k) = k_l*(r_alpha+r_cx/2.0D0)+k_u*(r_alpha-r_cy/2.0D0)
                        async_d(n_per_px,n_per_py,k) = k_r*(r_alpha-r_cx/2.0D0)+k_u*(r_alpha-r_cy/2.0D0)
                enddo

                ! Faces (6 planes)
                do i = 2,n_per_px-1
                    do j = 2, n_per_py-1
                        async_d(i,j,1) = k_b*(r_alpha+r_cz/2.0D0)
                        async_d(i,j,n_per_pz) = k_f*(r_alpha-r_cz/2.0D0)
                    enddo
                    do k = 2, n_per_pz-1
                        async_d(i,1,k) = k_d*(r_alpha+r_cy/2.0D0)
                        async_d(i,n_per_py,k) = k_u*(r_alpha-r_cy/2.0D0)
                    enddo
                enddo
                do j = 2,n_per_py-1
                    do k = 2, n_per_pz-1
                        async_d(1,j,k) = k_l*(r_alpha+r_cx/2.0D0)
                        async_d(n_per_px,j,k) = k_r*(r_alpha-r_cx/2.0D0)
                    enddo
                enddo

                do k = 1, n_per_pz
                    do j = 1, n_per_py
                        do i = 1, n_per_px
                            alpha_new = alpha/(1.0D0-async_d(i,j,k))
                            cx_new = cx/(1.0D0-async_d(i,j,k))
                            cy_new = cy/(1.0D0-async_d(i,j,k))
                            cz_new = cz/(1.0D0-async_d(i,j,k))

                            ! Proxy-Equation (Corrected Equation) Solution
                            call Sync_Scheme_3D(u_aproc(i,j,k,mod(t_par+1,L+1)),u_aproc(i-1,j,k,mod(t_par,L+1)),u_aproc(i,j-1,k,mod(t_par,L+1)), &
                                                u_aproc(i,j,k-1,mod(t_par,L+1)),u_aproc(i,j,k,mod(t_par,L+1)),u_aproc(i+1,j,k,mod(t_par,L+1)), &
                                                u_aproc(i,j+1,k,mod(t_par,L+1)),u_aproc(i,j,k+1,mod(t_par,L+1)),alpha_new,cx_new,cy_new,cz_new, &
                                                dt,dx,dy,dz)
                            u_aerr(i,j,k,mod(t_par+1,L+1)) = u_aproc(i,j,k,mod(t_par+1,L+1)) - u_exact(i,j,k)

                            ! Non-Corrected Equation Solution
                            call Sync_Scheme_3D(u_proc(i,j,k,mod(t_par+1,L+1)),u_proc(i-1,j,k,mod(t_par,L+1)),u_proc(i,j-1,k,mod(t_par,L+1)), &
                                                u_proc(i,j,k-1,mod(t_par,L+1)),u_proc(i,j,k,mod(t_par,L+1)),u_proc(i+1,j,k,mod(t_par,L+1)), &
                                                u_proc(i,j+1,k,mod(t_par,L+1)),u_proc(i,j,k+1,mod(t_par,L+1)),alpha,cx,cy,cz, &
                                                dt,dx,dy,dz)
                            u_err(i,j,k,mod(t_par+1,L+1)) = u_proc(i,j,k,mod(t_par+1,L+1)) - u_exact(i,j,k)
                        enddo
                    enddo
                enddo
        enddo
                if(taskid == 0) then
                    write(*,*) 'time loop ends'
                endif
                do k = 1, n_per_pz
                    do j = 1, n_per_py
                        do i = 1, n_per_px
                            Avg_aerr(k_par,n) = Avg_aerr(k_par,n) + dabs(u_aerr(i,j,k,mod(T,L+1)))
                            Avg_err(k_par,n) = Avg_err(k_par,n) + dabs(u_err(i,j,k,mod(T,L+1)))
                        enddo
                    enddo
                enddo
        enddo
                if(taskid == 0) then
                    write(*,*) 'angle loop ends'
                endif
                Avg_aerr(k_par,n) = Avg_aerr(k_par,n)/(nx(n)*ny(n)*nz(n)*angle_num)
                Avg_err(k_par,n) = Avg_err(k_par,n)/(nx(n)*ny(n)*nz(n)*angle_num)
                ref2(n) = (1.0D0/nx(n))*(1.0D0/nx(n))
                ref1(n) = (1.0D0/nx(n))
                if(taskid == 0) then
                        write(*,*) 'ref1 = ', ref1, 'ref2 = ', ref2
                endif
        enddo
                if(taskid == 0) then
                    write(*,*) 'kappa loop ends'
                endif
                deallocate(u_exact)
                deallocate(x)
                deallocate(y)
                deallocate(z)
                deallocate(u_proc)
                deallocate(u_err)
                deallocate(u_aproc)
                deallocate(u_aerr)
                deallocate(async_d)

                Call MPI_TYPE_FREE(back_front, error)
                Call MPI_TYPE_FREE(left_right, error)
                Call MPI_TYPE_FREE(down_up, error)
                Call MPI_TYPE_FREE(temp1, error)

                dt = 0.0D0
                dx = 0.0D0
                dy = 0.0D0
                dz = 0.0D0
        enddo

        call MPI_REDUCE(Avg_aerr, Avg_all_aerr, k_num*n_num, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, error)
        call MPI_REDUCE(Avg_err, Avg_all_err, k_num*n_num, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, error)

        if(taskid == 0) then
                ! Write Order Data (Corrected)
                write (filename, '( "Order_Async_corr_prob" , I2.2, ".dat" )' )  prob_par
                OPEN(UNIT = 550+prob_par, FILE = fileplace//filename, ACTION="write")
                do k_par = 1,K_num
                    do n = 1, n_num
                        write(550+prob_par,'(I1,X,I4,X,3(X,f21.14))') K_x, nx(n), Avg_all_aerr(k_par,n), ref1(n), ref2(n)
                    enddo
                enddo

                ! Write Order Data (Non-Corrected)
                write (filename, '( "Order_Async_non_corr_prob" , I2.2, ".dat" )' )  prob_par
                OPEN(UNIT = 600+prob_par, FILE = fileplace//filename, ACTION="write")
                do k_par = 1,K_num
                    do n = 1,n_num
                        write(600+prob_par,'(I1,X,I4,X,3(X,f21.14))') K_x, nx(n), Avg_all_err(k_par,n), ref1(n), ref2(n)
                    enddo
                enddo

                close(unit = 550+prob_par)
                close(unit = 600+prob_par)
        endif
        enddo

        call MPI_BARRIER(MPI_COMM_WORLD, error)
        write(*,*) 'Normal End of Execution'
        call MPI_FINALIZE(error)
end program
