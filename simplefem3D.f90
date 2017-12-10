!==============================================!
!                                              !
! C. thieulot ; December 2017                  !
!                                              !
!==============================================!
                                               !
program simpleFEM_2D                           !
                                               !
implicit none                                  !
                                               !
include 'mpif.h'                               !
include 'dmumps_struc.h'                       !
                                               !
integer, parameter :: m=8                      ! number of nodes which constitute an element
integer, parameter :: ndof=3                   ! number of dofs per node
integer nnx                                    ! number of grid points in the x direction
integer nny                                    ! number of grid points in the y direction
integer nnz                                    ! number of grid points in the y direction
integer np                                     ! number of grid points
integer nelx                                   ! number of elements in the x direction
integer nely                                   ! number of elements in the y direction
integer nelz                                   ! number of elements in the y direction
integer nel                                    ! number of elements
integer Nfem                                   ! size of the FEM matrix 
integer, dimension(:,:), allocatable :: icon   ! connectivity array
                                               !
integer i1,i2,i3,i,j,k,iel,counter,iq,jq,kq    !
integer ik,ikk,jkk,m1,k1,k2,ierr               ! integer parameters for loops
integer inode,idof,ic,iii,LELTVAR,NA_ELT       ! and bookkeeping
integer counter_mumps,ii,ij,iproc,nproc        !
                                               !  
real(8) Lx,Ly,Lz                               ! size of the numerical domain
real(8) viscosity                              ! dynamic viscosity $\mu$ of the material
real(8) density                                ! mass density $\rho$ of the material
real(8) gx,gy,gz                               ! gravity acceleration
real(8) penalty                                ! penalty parameter lambda
real(8), dimension(:), allocatable :: x,y,z    ! node coordinates arrays
real(8), dimension(:), allocatable :: u,v,w    ! node velocity arrays
real(8), dimension(:), allocatable :: press    ! pressure 
real(8), dimension(:), allocatable :: rho      ! density 
real(8), dimension(:), allocatable :: mu       ! viscosity
real(8), dimension(:), allocatable :: bc_val   ! array containing bc values
                                               !
real(8) rq,sq,tq,weightq                       ! local coordinate and weight of qpoint
real(8) xq,yq,zq                               ! global coordinate of qpoint
real(8) uq,vq,wq                               ! velocity at qpoint
real(8) exxq,eyyq,ezzq,exyq,exzq,eyzq          ! strain-rate components at qpoint  
real(8) Ael(m*ndof,m*ndof)                     ! elemental FEM matrix
real(8) Bel(m*ndof)                            ! elemental right hand side
real(8) N(m),dNdx(m),dNdy(m),dNdz(m)           !
real(8) dNdr(m),dNds(m),dNdt(m)                ! shape fcts and derivatives
real(8) jcob                                   ! determinant of jacobian matrix
real(8) jcb(3,3)                               ! jacobian matrix
real(8) jcbi(3,3)                              ! inverse of jacobian matrix
real(8) Bmat(6,ndof*m)                         ! B matrix
real(8), dimension(6,6) :: Kmat                ! K matrix 
real(8), dimension(6,6) :: Cmat                ! C matrix
real(8) fixt, Aref                             ! variables for applying b.c.
real(8), parameter :: epsiloon=1.d-10          ! tiny number
real(8) t1,t2,t3,t4                            !
                                               !
logical, dimension(:), allocatable :: bc_fix   ! prescribed b.c. array
                                               !
type(dmumps_struc) idV                         ! MUMPS data structure
                                               !
real(8), external :: rhofct,mufct              !
                                               !
!==============================================!
                                               !
CALL mpi_init(ierr)                            !  
call mpi_comm_size (mpi_comm_world,nproc,ierr) !
call mpi_comm_rank (mpi_comm_world,iproc,ierr) !
                                               !
!==============================================!
! initialise MUMPS                             !
!==============================================!
                                               !
idV%COMM = MPI_COMM_WORLD                      ! Define a communicator for the package 
idV%SYM = 1                                    ! Ask for symmetric matrix storage 
idV%par=1                                      ! Host working 
idV%JOB = -1                                   ! Initialize an instance of the package 
call DMUMPS(idV)                               ! MUMPS initialisation
                                               !
!==============================================!
!=====[setup]==================================!
!==============================================!
! grid counts np=nnx*nny nodes                 !
! and nel=nelx*nely elements.                  !
!==============================================!

if (iproc==0) open (unit=1020,file='OUT/mumps_output_info')

do nnx=4,50,2 ! loop over resolutions

   nny=nnx
   nnz=nnx

   Lx=1.d0
   Ly=1.d0
   Lz=1.d0

   np=nnx*nny*nnz

   nelx=nnx-1
   nely=nny-1
   nelz=nnz-1

   nel=nelx*nely*nelz

   penalty=1.d7

   Nfem=np*ndof

   Kmat=0.d0
   Kmat(1,1)=1.d0 ; Kmat(1,2)=1.d0 ; Kmat(1,3)=1.d0  
   Kmat(2,1)=1.d0 ; Kmat(2,2)=1.d0 ; Kmat(2,3)=1.d0  
   Kmat(3,1)=1.d0 ; Kmat(3,2)=1.d0 ; Kmat(3,3)=1.d0  

   Cmat=0.d0
   Cmat(1,1)=2.d0  ; Cmat(4,4)=1.d0
   Cmat(2,2)=2.d0  ; Cmat(5,5)=1.d0
   Cmat(3,3)=2.d0  ; Cmat(6,6)=1.d0

   gx=0.d0
   gy=0.d0
   gz=-1.d0

   write(*,'(a)') '======================================================================'
   write(*,'(a,3i5,a)') '=====[grid is ',nnx,nny,nnz,']=========================================='
   write(*,'(a)') '======================================================================'

   !==============================================!
   !===[allocate memory]==========================!
   !==============================================!
   call cpu_time(t3)

   allocate(x(np))
   allocate(y(np))
   allocate(z(np))
   allocate(u(np))
   allocate(v(np))
   allocate(w(np))
   allocate(icon(m,nel))
   allocate(bc_fix(Nfem))
   allocate(bc_val(Nfem))
   allocate(press(nel))
   allocate(rho(nel))
   allocate(mu(nel))

   call cpu_time(t4) ; write(*,*) Nfem,'allocate time      ',t4-t3

   !==============================================!
   !===[grid points setup]========================!
   !==============================================!

   call cpu_time(t3)

   counter=0
   do i=0,nelx
   do j=0,nely
   do k=0,nelz
      counter=counter+1
      x(counter)=dble(i)*Lx/dble(nelx)
      y(counter)=dble(j)*Ly/dble(nely)
      z(counter)=dble(k)*Lz/dble(nelz)
   end do
   end do
   end do

   call cpu_time(t4) ; write(*,*) Nfem,'grid setup time    ',t4-t3

   !==============================================!
   !===[connectivity]=============================!
   !==============================================!
   ! icon(1:8,iel) contains the node numbers of 
   ! element iel

   call cpu_time(t3)

   counter=0
   do i=1,nelx
   do j=1,nely
   do k=1,nelz
      counter=counter+1
      icon(1,counter)=nny*nnz*(i-1)+nnz*(j-1)+k
      icon(2,counter)=nny*nnz*(i  )+nnz*(j-1)+k
      icon(3,counter)=nny*nnz*(i  )+nnz*(j  )+k
      icon(4,counter)=nny*nnz*(i-1)+nnz*(j  )+k
      icon(5,counter)=nny*nnz*(i-1)+nnz*(j-1)+k+1
      icon(6,counter)=nny*nnz*(i  )+nnz*(j-1)+k+1
      icon(7,counter)=nny*nnz*(i  )+nnz*(j  )+k+1
      icon(8,counter)=nny*nnz*(i-1)+nnz*(j  )+k+1
   end do
   end do
   end do

   call cpu_time(t4) ; write(*,*) Nfem,'grid icon time     ',t4-t3

   !==============================================!
   ! MUMPS arrays
   !==============================================!

   call cpu_time(t3)

   idV%N=Nfem                                     ! total number of degrees of freedom, size of FEM matrix
   idV%NELT=nel                                   ! number of elements
   LELTVAR=nel*(m*ndof)                           ! nb of elts X size of elemental matrix
   NA_ELT=nel*(m*ndof)*(m*ndof+1)/2               ! nb of elts X nb of nbs in elemental matrix !NEW

   allocate(idV%A_ELT (NA_ELT))   
   allocate(idV%RHS   (idV%N))    

   if (iproc==0) then

      allocate(idV%ELTPTR(idV%NELT+1)) 
      allocate(idV%ELTVAR(LELTVAR))     
                     
      do i=1,nel                                  !
         idV%ELTPTR(i)=1+(i-1)*(ndof*m)           ! building ELTPTR array
      end do                                      !
      idV%ELTPTR(i)=1+nel*(ndof*m)                !
 
      counter=0                                   !
      do ic=1,nel                                 !
         do k=1,m                                 !
            inode=icon(k,ic)                      !
            do idof=1,ndof                        !
               iii=(inode-1)*ndof+idof            !
               counter=counter+1                  !
               idV%ELTVAR(counter)=iii            ! building ELTVAR
            end do                                !
         end do                                   !
      end do                                      !

   end if ! iproc

   idV%ICNTL(3) = 1020                            !
   idV%ICNTL(4) = 4                               !

   call cpu_time(t4) ; write(*,*) Nfem,'mumps arrays time  ',t4-t3

   !==============================================!
   !=====[define bc]==============================!
   !==============================================!

   call cpu_time(t3)

   bc_fix=.false.

   do i=1,np
      if (x(i).lt.epsiloon) then
         bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
         !bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
         !bc_fix((i-1)*ndof+3)=.true. ; bc_val((i-1)*ndof+3)=0.d0
      endif
      if (x(i).gt.(Lx-epsiloon)) then
         bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
         !bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
         !bc_fix((i-1)*ndof+3)=.true. ; bc_val((i-1)*ndof+3)=0.d0
      endif
      if (y(i).lt.epsiloon) then
         !bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
         bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
         !bc_fix((i-1)*ndof+3)=.true. ; bc_val((i-1)*ndof+3)=0.d0
      endif
      if (y(i).gt.(Ly-epsiloon) ) then
         !bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
         bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
         !bc_fix((i-1)*ndof+3)=.true. ; bc_val((i-1)*ndof+3)=0.d0
      endif
      if (z(i).lt.epsiloon) then
         !bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
         !bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
         bc_fix((i-1)*ndof+3)=.true. ; bc_val((i-1)*ndof+3)=0.d0
      endif
      if (z(i).gt.(Lz-epsiloon) ) then
         !bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
         !bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
         bc_fix((i-1)*ndof+3)=.true. ; bc_val((i-1)*ndof+3)=0.d0
      endif
   end do

   call cpu_time(t4) ; write(*,*) Nfem,'bc time            ',t4-t3

   !==============================================!
   !=====[build FE matrix]========================!
   !==============================================!

   call cpu_time(t3)

   idV%RHS=0.d0
   idV%A_ELT=0.d0
   counter_mumps=0

   do iel=1,nel

      Ael=0.d0
      Bel=0.d0

      do iq=-1,1,2
      do jq=-1,1,2
      do kq=-1,1,2

         rq=iq/sqrt(3.d0)
         sq=jq/sqrt(3.d0)
         tq=kq/sqrt(3.d0)
         weightq=1.d0*1.d0*1.d0

         N(1)=0.125d0*(1.d0-rq)*(1.d0-sq)*(1.d0-tq)
         N(2)=0.125d0*(1.d0+rq)*(1.d0-sq)*(1.d0-tq)
         N(3)=0.125d0*(1.d0+rq)*(1.d0+sq)*(1.d0-tq)
         N(4)=0.125d0*(1.d0-rq)*(1.d0+sq)*(1.d0-tq)
         N(5)=0.125d0*(1.d0-rq)*(1.d0-sq)*(1.d0+tq)
         N(6)=0.125d0*(1.d0+rq)*(1.d0-sq)*(1.d0+tq)
         N(7)=0.125d0*(1.d0+rq)*(1.d0+sq)*(1.d0+tq)
         N(8)=0.125d0*(1.d0-rq)*(1.d0+sq)*(1.d0+tq)

         dNdr(1)= - 0.125d0*(1.d0-sq)*(1.d0-tq)
         dNdr(2)= + 0.125d0*(1.d0-sq)*(1.d0-tq)
         dNdr(3)= + 0.125d0*(1.d0+sq)*(1.d0-tq)
         dNdr(4)= - 0.125d0*(1.d0+sq)*(1.d0-tq)
         dNdr(5)= - 0.125d0*(1.d0-sq)*(1.d0+tq)
         dNdr(6)= + 0.125d0*(1.d0-sq)*(1.d0+tq)
         dNdr(7)= + 0.125d0*(1.d0+sq)*(1.d0+tq)
         dNdr(8)= - 0.125d0*(1.d0+sq)*(1.d0+tq)

         dNds(1)= - 0.125d0*(1.d0-rq)*(1.d0-tq)
         dNds(2)= - 0.125d0*(1.d0+rq)*(1.d0-tq)
         dNds(3)= + 0.125d0*(1.d0+rq)*(1.d0-tq)
         dNds(4)= + 0.125d0*(1.d0-rq)*(1.d0-tq)
         dNds(5)= - 0.125d0*(1.d0-rq)*(1.d0+tq)
         dNds(6)= - 0.125d0*(1.d0+rq)*(1.d0+tq)
         dNds(7)= + 0.125d0*(1.d0+rq)*(1.d0+tq)
         dNds(8)= + 0.125d0*(1.d0-rq)*(1.d0+tq)

         dNdt(1)= - 0.125d0*(1.d0-rq)*(1.d0-sq)
         dNdt(2)= - 0.125d0*(1.d0+rq)*(1.d0-sq)
         dNdt(3)= - 0.125d0*(1.d0+rq)*(1.d0+sq)
         dNdt(4)= - 0.125d0*(1.d0-rq)*(1.d0+sq)
         dNdt(5)= + 0.125d0*(1.d0-rq)*(1.d0-sq)
         dNdt(6)= + 0.125d0*(1.d0+rq)*(1.d0-sq)
         dNdt(7)= + 0.125d0*(1.d0+rq)*(1.d0+sq)
         dNdt(8)= + 0.125d0*(1.d0-rq)*(1.d0+sq)

         jcb=0.d0
         do k=1,m
            jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
            jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
            jcb(1,3)=jcb(1,3)+dNdr(k)*z(icon(k,iel))
            jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
            jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
            jcb(2,3)=jcb(2,3)+dNds(k)*z(icon(k,iel))
            jcb(3,1)=jcb(3,1)+dNdt(k)*x(icon(k,iel))
            jcb(3,2)=jcb(3,2)+dNdt(k)*y(icon(k,iel))
            jcb(3,3)=jcb(3,3)+dNdt(k)*z(icon(k,iel))
         enddo

         jcob=jcb(1,1)*jcb(2,2)*jcb(3,3) &
             +jcb(1,2)*jcb(2,3)*jcb(3,1) &
             +jcb(2,1)*jcb(3,2)*jcb(1,3) &
             -jcb(1,3)*jcb(2,2)*jcb(3,1) &
             -jcb(1,2)*jcb(2,1)*jcb(3,3) &
             -jcb(2,3)*jcb(3,2)*jcb(1,1)

         jcbi(1,1)=(jcb(2,2)*jcb(3,3)-jcb(2,3)*jcb(3,2))/jcob
         jcbi(2,1)=(jcb(2,3)*jcb(3,1)-jcb(2,1)*jcb(3,3))/jcob
         jcbi(3,1)=(jcb(2,1)*jcb(3,2)-jcb(2,2)*jcb(3,1))/jcob
         jcbi(1,2)=(jcb(1,3)*jcb(3,2)-jcb(1,2)*jcb(3,3))/jcob
         jcbi(2,2)=(jcb(1,1)*jcb(3,3)-jcb(1,3)*jcb(3,1))/jcob
         jcbi(3,2)=(jcb(1,2)*jcb(3,1)-jcb(1,1)*jcb(3,2))/jcob
         jcbi(1,3)=(jcb(1,2)*jcb(2,3)-jcb(1,3)*jcb(2,2))/jcob
         jcbi(2,3)=(jcb(1,3)*jcb(2,1)-jcb(1,1)*jcb(2,3))/jcob
         jcbi(3,3)=(jcb(1,1)*jcb(2,2)-jcb(1,2)*jcb(2,1))/jcob

         xq=0.d0 ; yq=0.d0 ; zq=0.d0
         uq=0.d0 ; vq=0.d0 ; wq=0.d0
         exxq=0.d0 ; eyyq=0.d0 ; ezzq=0.d0
         exyq=0.d0 ; exzq=0.d0 ; eyzq=0.d0
         do k=1,m
            xq=xq+N(k)*x(icon(k,iel))
            yq=yq+N(k)*y(icon(k,iel))
            zq=zq+N(k)*z(icon(k,iel))
            uq=uq+N(k)*u(icon(k,iel))
            vq=vq+N(k)*v(icon(k,iel))
            wq=wq+N(k)*w(icon(k,iel))
            dNdx(k)=jcbi(1,1)*dNdr(k)&
                   +jcbi(1,2)*dNds(k)&
                   +jcbi(1,3)*dNdt(k)
            dNdy(k)=jcbi(2,1)*dNdr(k)&
                   +jcbi(2,2)*dNds(k)&
                   +jcbi(2,3)*dNdt(k)
            dNdz(k)=jcbi(3,1)*dNdr(k)&
                   +jcbi(3,2)*dNds(k)&
                   +jcbi(3,3)*dNdt(k)
            exxq=exxq+ dNdx(k)*u(icon(k,iel))
            eyyq=eyyq+ dNdy(k)*v(icon(k,iel))
            ezzq=ezzq+ dNdz(k)*w(icon(k,iel))
            exyq=exyq+ dNdx(k)*v(icon(k,iel)) *0.5d0 &
                     + dNdy(k)*u(icon(k,iel)) *0.5d0
            exzq=exzq+ dNdx(k)*w(icon(k,iel)) *0.5d0 &
                     + dNdz(k)*u(icon(k,iel)) *0.5d0
            eyzq=eyzq+ dNdy(k)*w(icon(k,iel)) *0.5d0 &
                     + dNdz(k)*v(icon(k,iel)) *0.5d0
         end do

         Bmat=0.d0
         do i=1,m
         i1=ndof*i-2
         i2=ndof*i-1
         i3=ndof*i
         Bmat(1,i1)=dNdx(i)
         Bmat(2,i2)=dNdy(i)
         Bmat(3,i3)=dNdz(i)
         Bmat(4,i1)=dNdy(i) ; Bmat(4,i2)=dNdx(i)
         Bmat(5,i1)=dNdz(i) ; Bmat(5,i3)=dNdx(i)
         Bmat(6,i2)=dNdz(i) ; Bmat(6,i3)=dNdy(i)
         end do

         viscosity=mufct(xq,yq,zq)
         Ael=Ael + matmul(transpose(Bmat),matmul(viscosity*Cmat,Bmat))*weightq*jcob

         do i=1,m
         i1=ndof*i-2
         i2=ndof*i-1
         i3=ndof*i
         density=rhofct(xq,yq,zq)
         Bel(i1)=Bel(i1)+N(i)*jcob*weightq*density*gx
         Bel(i2)=Bel(i2)+N(i)*jcob*weightq*density*gy
         Bel(i3)=Bel(i3)+N(i)*jcob*weightq*density*gz
         end do

      end do
      end do
      end do

      ! 1 point integration

      rq=0.d0
      sq=0.d0
      tq=0.d0
      weightq=2.d0*2.d0*2.d0

      N(1)=0.125d0*(1.d0-rq)*(1.d0-sq)*(1.d0-tq)
      N(2)=0.125d0*(1.d0+rq)*(1.d0-sq)*(1.d0-tq)
      N(3)=0.125d0*(1.d0+rq)*(1.d0+sq)*(1.d0-tq)
      N(4)=0.125d0*(1.d0-rq)*(1.d0+sq)*(1.d0-tq)
      N(5)=0.125d0*(1.d0-rq)*(1.d0-sq)*(1.d0+tq)
      N(6)=0.125d0*(1.d0+rq)*(1.d0-sq)*(1.d0+tq)
      N(7)=0.125d0*(1.d0+rq)*(1.d0+sq)*(1.d0+tq)
      N(8)=0.125d0*(1.d0-rq)*(1.d0+sq)*(1.d0+tq)

      dNdr(1)= - 0.125d0*(1.d0-sq)*(1.d0-tq)
      dNdr(2)= + 0.125d0*(1.d0-sq)*(1.d0-tq)
      dNdr(3)= + 0.125d0*(1.d0+sq)*(1.d0-tq)
      dNdr(4)= - 0.125d0*(1.d0+sq)*(1.d0-tq)
      dNdr(5)= - 0.125d0*(1.d0-sq)*(1.d0+tq)
      dNdr(6)= + 0.125d0*(1.d0-sq)*(1.d0+tq)
      dNdr(7)= + 0.125d0*(1.d0+sq)*(1.d0+tq)
      dNdr(8)= - 0.125d0*(1.d0+sq)*(1.d0+tq)

      dNds(1)= - 0.125d0*(1.d0-rq)*(1.d0-tq)
      dNds(2)= - 0.125d0*(1.d0+rq)*(1.d0-tq)
      dNds(3)= + 0.125d0*(1.d0+rq)*(1.d0-tq)
      dNds(4)= + 0.125d0*(1.d0-rq)*(1.d0-tq)
      dNds(5)= - 0.125d0*(1.d0-rq)*(1.d0+tq)
      dNds(6)= - 0.125d0*(1.d0+rq)*(1.d0+tq)
      dNds(7)= + 0.125d0*(1.d0+rq)*(1.d0+tq)
      dNds(8)= + 0.125d0*(1.d0-rq)*(1.d0+tq)

      dNdt(1)= - 0.125d0*(1.d0-rq)*(1.d0-sq)
      dNdt(2)= - 0.125d0*(1.d0+rq)*(1.d0-sq)
      dNdt(3)= - 0.125d0*(1.d0+rq)*(1.d0+sq)
      dNdt(4)= - 0.125d0*(1.d0-rq)*(1.d0+sq)
      dNdt(5)= + 0.125d0*(1.d0-rq)*(1.d0-sq)
      dNdt(6)= + 0.125d0*(1.d0+rq)*(1.d0-sq)
      dNdt(7)= + 0.125d0*(1.d0+rq)*(1.d0+sq)
      dNdt(8)= + 0.125d0*(1.d0-rq)*(1.d0+sq)

      jcb=0.d0
      do k=1,8
      jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
      jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
      jcb(1,3)=jcb(1,3)+dNdr(k)*z(icon(k,iel))
      jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
      jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      jcb(2,3)=jcb(2,3)+dNds(k)*z(icon(k,iel))
      jcb(3,1)=jcb(3,1)+dNdt(k)*x(icon(k,iel))
      jcb(3,2)=jcb(3,2)+dNdt(k)*y(icon(k,iel))
      jcb(3,3)=jcb(3,3)+dNdt(k)*z(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)*jcb(3,3) &
          +jcb(1,2)*jcb(2,3)*jcb(3,1) &
          +jcb(2,1)*jcb(3,2)*jcb(1,3) &
          -jcb(1,3)*jcb(2,2)*jcb(3,1) &
          -jcb(1,2)*jcb(2,1)*jcb(3,3) &
          -jcb(2,3)*jcb(3,2)*jcb(1,1)

      jcbi(1,1)=(jcb(2,2)*jcb(3,3)-jcb(2,3)*jcb(3,2))/jcob
      jcbi(2,1)=(jcb(2,3)*jcb(3,1)-jcb(2,1)*jcb(3,3))/jcob
      jcbi(3,1)=(jcb(2,1)*jcb(3,2)-jcb(2,2)*jcb(3,1))/jcob
      jcbi(1,2)=(jcb(1,3)*jcb(3,2)-jcb(1,2)*jcb(3,3))/jcob
      jcbi(2,2)=(jcb(1,1)*jcb(3,3)-jcb(1,3)*jcb(3,1))/jcob
      jcbi(3,2)=(jcb(1,2)*jcb(3,1)-jcb(1,1)*jcb(3,2))/jcob
      jcbi(1,3)=(jcb(1,2)*jcb(2,3)-jcb(1,3)*jcb(2,2))/jcob
      jcbi(2,3)=(jcb(1,3)*jcb(2,1)-jcb(1,1)*jcb(2,3))/jcob
      jcbi(3,3)=(jcb(1,1)*jcb(2,2)-jcb(1,2)*jcb(2,1))/jcob

      do k=1,m
         dNdx(k)=jcbi(1,1)*dNdr(k)&
                +jcbi(1,2)*dNds(k)&
                +jcbi(1,3)*dNdt(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)&
                +jcbi(2,2)*dNds(k)&
                +jcbi(2,3)*dNdt(k)
         dNdz(k)=jcbi(3,1)*dNdr(k)&
                +jcbi(3,2)*dNds(k)&
                +jcbi(3,3)*dNdt(k)
      end do

      Bmat=0.d0
      do i=1,m
         i1=ndof*i-2
         i2=ndof*i-1
         i3=ndof*i
         Bmat(1,i1)=dNdx(i)
         Bmat(2,i2)=dNdy(i)
         Bmat(3,i3)=dNdz(i)
         Bmat(4,i1)=dNdy(i) ; Bmat(4,i2)=dNdx(i)
         Bmat(5,i1)=dNdz(i) ; Bmat(5,i3)=dNdx(i)
         Bmat(6,i2)=dNdz(i) ; Bmat(6,i3)=dNdy(i)
      end do

      Ael=Ael + matmul(transpose(Bmat),matmul(penalty*Kmat,Bmat))*weightq*jcob

      !======================================
      !=====[impose boundary conditions]=====
      !======================================
                
      do ii=1,m  
         inode=icon(ii,iel)  
         do k=1,ndof       
            ij=(inode-1)*ndof+k    
            if (bc_fix(ij)) then  
            fixt=bc_val(ij) 
            i=(ii-1)*ndof+k               
            Aref=Ael(i,i)                
            do j=1,m*ndof               
               Bel(j)=Bel(j)-Ael(j,i)*fixt 
               Ael(i,j)=0.d0              
               Ael(j,i)=0.d0             
            enddo                       
            Ael(i,i)=Aref              
            Bel(i)=Aref*fixt          
            endif    
         enddo      
      enddo        

      !=====================
      !=====[assemble]======
      !=====================

      do k1=1,m   
         ik=icon(k1,iel) 
         do i1=1,ndof  
            ikk=ndof*(k1-1)+i1    
            m1=ndof*(ik-1)+i1    
            do k2=1,m         
               do i2=1,ndof    
                  jkk=ndof*(k2-1)+i2  
                  if (jkk>=ikk) then 
                     counter_mumps=counter_mumps+1   
                     idV%A_ELT(counter_mumps)=Ael(ikk,jkk) 
                  end if                                  
               end do                                    
            end do                                      
            idV%RHS(m1)=idV%RHS(m1)+Bel(ikk)           
         end do                                  
      end do      

   end do

   call cpu_time(t4) ; write(*,*) Nfem,'make matrix        ',t4-t3

   !==============================================!
   !=====[solve system]===========================!
   !==============================================!

   idV%ICNTL(5) = 1                               ! elemental format
   idV%ICNTL(18) = 0                              ! the input matrix is centralized on the host

   call cpu_time(t3)

   call cpu_time(t1)
   idV%JOB = 1                                    ! analysis phase 
   CALL DMUMPS(idV)
   call cpu_time(t2) 
   write(*,*) Nfem,'analysis time ',t2-t1,'s'

   call cpu_time(t1)
   idV%JOB = 2                                    ! factorisation phase 
   CALL DMUMPS(idV)
   call cpu_time(t2)
   write(*,*) Nfem,'factor.  time ',t2-t1,'s'

   call cpu_time(t1)
   idV%JOB = 3                                    ! solve phase
   CALL DMUMPS(idV)
   call cpu_time(t2)
   write(*,*) Nfem,'solution time ',t2-t1,'s'

   call cpu_time(t4) ; write(*,*) Nfem,'total time    ',t4-t3,'s'

   write(*,*) Nfem,'estimated RAM needed for the factorization',idV%info(15),'Mb'
   write(*,*) Nfem,'effective RAM needed for the factorization',idV%info(22),'Mb'
   write(*,*) Nfem,'effective RAM needed for the solution     ',idV%info(26),'Mb'

   !if (nproc>1 .and. iproc==0) then
   !write(*,'(a,i7,a)') 'estimated RAM needed for the factorization:',idV%infog(17),'Mb (sum over all MPI processes)'
   !write(*,'(a,i7,a)') 'RAM allocated        for the factorization:',idV%infog(19),'Mb (sum over all MPI processes)'
   !write(*,'(a,i7,a)') 'effective RAM needed for the factorization:',idV%infog(22),'Mb (sum over all MPI processes)'
   !write(*,'(a,i7,a)') 'effective RAM needed for the solution     :',idV%infog(31),'Mb (sum over all MPI processes)'
   !end if

   !==============================================!
   !=====[transfer solution]======================!
   !==============================================!

   do i=1,np
      u(i)=idV%RHS((i-1)*ndof+1)
      v(i)=idV%RHS((i-1)*ndof+2)
      w(i)=idV%RHS((i-1)*ndof+3)
   end do

   !==============================================!
   !=====[retrieve pressure]======================!
   !==============================================!

   do iel=1,nel

      rq=0.d0
      sq=0.d0
      tq=0.d0

      N(1)=0.125d0*(1.d0-rq)*(1.d0-sq)*(1.d0-tq)
      N(2)=0.125d0*(1.d0+rq)*(1.d0-sq)*(1.d0-tq)
      N(3)=0.125d0*(1.d0+rq)*(1.d0+sq)*(1.d0-tq)
      N(4)=0.125d0*(1.d0-rq)*(1.d0+sq)*(1.d0-tq)
      N(5)=0.125d0*(1.d0-rq)*(1.d0-sq)*(1.d0+tq)
      N(6)=0.125d0*(1.d0+rq)*(1.d0-sq)*(1.d0+tq)
      N(7)=0.125d0*(1.d0+rq)*(1.d0+sq)*(1.d0+tq)
      N(8)=0.125d0*(1.d0-rq)*(1.d0+sq)*(1.d0+tq)

      dNdr(1)= - 0.125d0*(1.d0-sq)*(1.d0-tq)
      dNdr(2)= + 0.125d0*(1.d0-sq)*(1.d0-tq)
      dNdr(3)= + 0.125d0*(1.d0+sq)*(1.d0-tq)
      dNdr(4)= - 0.125d0*(1.d0+sq)*(1.d0-tq)
      dNdr(5)= - 0.125d0*(1.d0-sq)*(1.d0+tq)
      dNdr(6)= + 0.125d0*(1.d0-sq)*(1.d0+tq)
      dNdr(7)= + 0.125d0*(1.d0+sq)*(1.d0+tq)
      dNdr(8)= - 0.125d0*(1.d0+sq)*(1.d0+tq)

      dNds(1)= - 0.125d0*(1.d0-rq)*(1.d0-tq)
      dNds(2)= - 0.125d0*(1.d0+rq)*(1.d0-tq)
      dNds(3)= + 0.125d0*(1.d0+rq)*(1.d0-tq)
      dNds(4)= + 0.125d0*(1.d0-rq)*(1.d0-tq)
      dNds(5)= - 0.125d0*(1.d0-rq)*(1.d0+tq)
      dNds(6)= - 0.125d0*(1.d0+rq)*(1.d0+tq)
      dNds(7)= + 0.125d0*(1.d0+rq)*(1.d0+tq)
      dNds(8)= + 0.125d0*(1.d0-rq)*(1.d0+tq)

      dNdt(1)= - 0.125d0*(1.d0-rq)*(1.d0-sq)
      dNdt(2)= - 0.125d0*(1.d0+rq)*(1.d0-sq)
      dNdt(3)= - 0.125d0*(1.d0+rq)*(1.d0+sq)
      dNdt(4)= - 0.125d0*(1.d0-rq)*(1.d0+sq)
      dNdt(5)= + 0.125d0*(1.d0-rq)*(1.d0-sq)
      dNdt(6)= + 0.125d0*(1.d0+rq)*(1.d0-sq)
      dNdt(7)= + 0.125d0*(1.d0+rq)*(1.d0+sq)
      dNdt(8)= + 0.125d0*(1.d0-rq)*(1.d0+sq)

      jcb=0.d0
      do k=1,8
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(1,3)=jcb(1,3)+dNdr(k)*z(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
         jcb(2,3)=jcb(2,3)+dNds(k)*z(icon(k,iel))
         jcb(3,1)=jcb(3,1)+dNdt(k)*x(icon(k,iel))
         jcb(3,2)=jcb(3,2)+dNdt(k)*y(icon(k,iel))
         jcb(3,3)=jcb(3,3)+dNdt(k)*z(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)*jcb(3,3) &
          +jcb(1,2)*jcb(2,3)*jcb(3,1) &
          +jcb(2,1)*jcb(3,2)*jcb(1,3) &
          -jcb(1,3)*jcb(2,2)*jcb(3,1) &
          -jcb(1,2)*jcb(2,1)*jcb(3,3) &
          -jcb(2,3)*jcb(3,2)*jcb(1,1)

      jcbi(1,1)=(jcb(2,2)*jcb(3,3)-jcb(2,3)*jcb(3,2))/jcob
      jcbi(2,1)=(jcb(2,3)*jcb(3,1)-jcb(2,1)*jcb(3,3))/jcob
      jcbi(3,1)=(jcb(2,1)*jcb(3,2)-jcb(2,2)*jcb(3,1))/jcob
      jcbi(1,2)=(jcb(1,3)*jcb(3,2)-jcb(1,2)*jcb(3,3))/jcob
      jcbi(2,2)=(jcb(1,1)*jcb(3,3)-jcb(1,3)*jcb(3,1))/jcob
      jcbi(3,2)=(jcb(1,2)*jcb(3,1)-jcb(1,1)*jcb(3,2))/jcob
      jcbi(1,3)=(jcb(1,2)*jcb(2,3)-jcb(1,3)*jcb(2,2))/jcob
      jcbi(2,3)=(jcb(1,3)*jcb(2,1)-jcb(1,1)*jcb(2,3))/jcob
      jcbi(3,3)=(jcb(1,1)*jcb(2,2)-jcb(1,2)*jcb(2,1))/jcob

      do k=1,m
         dNdx(k)=jcbi(1,1)*dNdr(k)&
                +jcbi(1,2)*dNds(k)&
                +jcbi(1,3)*dNdt(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)&
                +jcbi(2,2)*dNds(k)&
                +jcbi(2,3)*dNdt(k)
         dNdz(k)=jcbi(3,1)*dNdr(k)&
                +jcbi(3,2)*dNds(k)&
                +jcbi(3,3)*dNdt(k)
      end do

      xq=0.d0
      yq=0.d0
      zq=0.d0
      exxq=0.d0
      eyyq=0.d0
      ezzq=0.d0
      do k=1,m
         xq=xq+N(k)*x(icon(k,iel))
         yq=yq+N(k)*y(icon(k,iel))
         zq=zq+N(k)*z(icon(k,iel))
         exxq=exxq+ dNdx(k)*u(icon(k,iel))
         eyyq=eyyq+ dNdy(k)*v(icon(k,iel))
         ezzq=ezzq+ dNdz(k)*w(icon(k,iel))
      end do

      press(iel)=-penalty*(exxq+eyyq+ezzq)

      rho(iel)=rhofct(xq,yq,zq)
      mu(iel)=mufct(xq,yq,zq)

   end do

   !======================
   !open(unit=123,file='OUT/solution_u.dat',status='replace')
   !open(unit=234,file='OUT/solution_v.dat',status='replace')
   !do i=1,np
   !   write(123,'(5f20.10)') x(i),y(i),u(i),uth(x(i),y(i)),u(i)-
   !   write(234,'(5f20.10)') x(i),y(i),v(i),vth(x(i),y(i)),
   !end do
   !close(123)
   !close(234)
   !======================

   if (iproc==0) call output_for_paraview3D (np,nel,x,y,z,u,v,w,press,icon,rho,mu)

   deallocate(x,y,z)
   deallocate(u,v,w)
   deallocate(icon)
   deallocate(bc_fix)
   deallocate(bc_val)
   deallocate(press)
   deallocate(rho)
   deallocate(mu)
   deallocate(idV%A_ELT)   
   deallocate(idV%RHS)    
   deallocate(idV%ELTPTR) 
   deallocate(idV%ELTVAR)     

end do ! nnx 


if (iproc==0) close (1020)

end program

!==============================================!
!==============================================!
!==============================================!

function rhofct (x,y,z)
implicit none
real(8) rhofct,x,y,z

if (sqrt((x-0.5)**2+(y-0.5)**2+(z-0.5)**2)< 1./8.) then
   rhofct=2.d0
else
   rhofct=1.d0
end if

end function

function mufct (x,y,z)
implicit none
real(8) mufct,x,y,z

if (sqrt((x-0.5)**2+(y-0.5)**2+(z-0.5)**2)< 1./8.) then
   mufct=100.d0
else
   mufct=1.d0
end if

end function









