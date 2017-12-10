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
integer, parameter :: m=4                      ! number of nodes which constitute an element
integer, parameter :: ndof=2                   ! number of dofs per node
integer nnx                                    ! number of grid points in the x direction
integer nny                                    ! number of grid points in the y direction
integer np                                     ! number of grid points
integer nelx                                   ! number of elements in the x direction
integer nely                                   ! number of elements in the y direction
integer nel                                    ! number of elements
integer Nfem                                   ! size of the FEM matrix 
integer, dimension(:,:), allocatable :: icon   ! connectivity array
                                               !
integer i1,i2,i,j,k,iel,counter,iq,jq,kkk      !
integer ik,ikk,jkk,m1,k1,k2,ierr               ! integer parameters for loops
integer inode,idof,ic,iii,LELTVAR,NA_ELT       ! and bookkeeping
integer counter_mumps,ii,ij,iproc,nproc        !
                                               !  
real(8) Lx,Ly                                  ! size of the numerical domain
real(8) viscosity                              ! dynamic viscosity $\mu$ of the material
real(8) density                                ! mass density $\rho$ of the material
real(8) gx,gy                                  ! gravity acceleration
real(8) penalty                                ! penalty parameter lambda
real(8), dimension(:), allocatable :: x,y      ! node coordinates arrays
real(8), dimension(:), allocatable :: u,v      ! node velocity arrays
real(8), dimension(:), allocatable :: press    ! pressure 
real(8), dimension(:), allocatable :: bc_val   ! array containing bc values
                                               !
real(8), external :: b1,b2,uth,vth,pth         ! body force and analytical solution
real(8) rq,sq,wq                               ! local coordinate and weight of qpoint
real(8) xq,yq                                  ! global coordinate of qpoint
real(8) uq,vq                                  ! velocity at qpoint
real(8) exxq,eyyq,exyq                         ! strain-rate components at qpoint  
real(8) Ael(m*ndof,m*ndof)                     ! elemental FEM matrix
real(8) Bel(m*ndof)                            ! elemental right hand side
real(8) N(m),dNdx(m),dNdy(m),dNdr(m),dNds(m)   ! shape fcts and derivatives
real(8) jcob                                   ! determinant of jacobian matrix
real(8) jcb(2,2)                               ! jacobian matrix
real(8) jcbi(2,2)                              ! inverse of jacobian matrix
real(8) Bmat(3,ndof*m)                         ! B matrix
real(8), dimension(3,3) :: Kmat                ! K matrix 
real(8), dimension(3,3) :: Cmat                ! C matrix
real(8) fixt, Aref                             ! variables for applying b.c.
real(8), parameter :: epsiloon=1.d-10          ! tiny number
real(8) t1,t2,t3,t4                            !
                                               !
logical, dimension(:), allocatable :: bc_fix   ! prescribed b.c. array
                                               !
type(dmumps_struc) idV                         ! MUMPS data structure
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

do kkk=4,7!22 ! loop over resolutions

   nnx=int(sqrt(2.d0**kkk))
   nny=nnx

   Lx=1.d0
   Ly=1.d0

   np=nnx*nny

   nelx=nnx-1
   nely=nny-1

   nel=nelx*nely

   penalty=1.d7

   viscosity=1.d0
   density=1.d0

   Nfem=np*ndof

   Kmat(1,1)=1.d0 ; Kmat(1,2)=1.d0 ; Kmat(1,3)=0.d0  
   Kmat(2,1)=1.d0 ; Kmat(2,2)=1.d0 ; Kmat(2,3)=0.d0  
   Kmat(3,1)=0.d0 ; Kmat(3,2)=0.d0 ; Kmat(3,3)=0.d0  

   Cmat(1,1)=2.d0 ; Cmat(1,2)=0.d0 ; Cmat(1,3)=0.d0  
   Cmat(2,1)=0.d0 ; Cmat(2,2)=2.d0 ; Cmat(2,3)=0.d0  
   Cmat(3,1)=0.d0 ; Cmat(3,2)=0.d0 ; Cmat(3,3)=1.d0  

   write(*,'(a)') '======================================================================'
   write(*,'(a,2i5,a)') '=====[grid is ',nnx,nny,']============================================='
   write(*,'(a)') '======================================================================'

   !==============================================!
   !===[allocate memory]==========================!
   !==============================================!
   call cpu_time(t3)

   allocate(x(np))
   allocate(y(np))
   allocate(u(np))
   allocate(v(np))
   allocate(icon(m,nel))
   allocate(bc_fix(Nfem))
   allocate(bc_val(Nfem))
   allocate(press(nel))

   call cpu_time(t4) ; write(*,*) Nfem,'allocate time      ',t4-t3

   !==============================================!
   !===[grid points setup]========================!
   !==============================================!

   call cpu_time(t3)

   counter=0
   do j=0,nely
   do i=0,nelx
      counter=counter+1
      x(counter)=dble(i)*Lx/dble(nelx)
      y(counter)=dble(j)*Ly/dble(nely)
   end do
   end do

   call cpu_time(t4) ; write(*,*) Nfem,'grid setup time    ',t4-t3

   !==============================================!
   !===[connectivity]=============================!
   !==============================================!
   ! icon(1:4,iel) contains the node numbers of 
   ! element iel

   call cpu_time(t3)

   counter=0
   do j=1,nely
   do i=1,nelx
      counter=counter+1
      icon(1,counter)=i+(j-1)*(nelx+1)
      icon(2,counter)=i+1+(j-1)*(nelx+1)
      icon(3,counter)=i+1+j*(nelx+1)
      icon(4,counter)=i+j*(nelx+1)
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
         bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
      endif
      if (x(i).gt.(Lx-epsiloon)) then
         bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
         bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
      endif
      if (y(i).lt.epsiloon) then
         bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
         bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
      endif
      if (y(i).gt.(Ly-epsiloon) ) then
         bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
         bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
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

         rq=iq/sqrt(3.d0)
         sq=jq/sqrt(3.d0)
         wq=1.d0*1.d0

         N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
         N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
         N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
         N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

         dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
         dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
         dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
         dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

         jcb=0.d0
         do k=1,m
            jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
            jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
            jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
            jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
         enddo

         jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

         jcbi(1,1)=    jcb(2,2) /jcob
         jcbi(1,2)=  - jcb(1,2) /jcob
         jcbi(2,1)=  - jcb(2,1) /jcob
         jcbi(2,2)=    jcb(1,1) /jcob

         xq=0.d0
         yq=0.d0
         uq=0.d0
         vq=0.d0
         exxq=0.d0
         eyyq=0.d0
         exyq=0.d0
         do k=1,m
            xq=xq+N(k)*x(icon(k,iel))
            yq=yq+N(k)*y(icon(k,iel))
            uq=uq+N(k)*u(icon(k,iel))
            vq=vq+N(k)*v(icon(k,iel))
            dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
            dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
            exxq=exxq+ dNdx(k)*u(icon(k,iel))
            eyyq=eyyq+ dNdy(k)*v(icon(k,iel))
            exyq=exyq+ dNdx(k)*v(icon(k,iel)) *0.5d0 &
                     + dNdy(k)*u(icon(k,iel)) *0.5d0
         end do

         do i=1,m
            i1=2*i-1
            i2=2*i
            Bmat(1,i1)=dNdx(i) ; Bmat(1,i2)=0.d0
            Bmat(2,i1)=0.d0    ; Bmat(2,i2)=dNdy(i)
            Bmat(3,i1)=dNdy(i) ; Bmat(3,i2)=dNdx(i)
         end do

         Ael=Ael + matmul(transpose(Bmat),matmul(viscosity*Cmat,Bmat))*wq*jcob

         do i=1,m
            i1=2*i-1
            i2=2*i
            !Bel(i1)=Bel(i1)+N(i)*jcob*wq*density*gx
            !Bel(i2)=Bel(i2)+N(i)*jcob*wq*density*gy
            Bel(i1)=Bel(i1)+N(i)*jcob*wq*b1(xq,yq)
            Bel(i2)=Bel(i2)+N(i)*jcob*wq*b2(xq,yq)
         end do

      end do
      end do

      ! 1 point integration

      rq=0.d0
      sq=0.d0
      wq=2.d0*2.d0

      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      jcbi(1,1)=    jcb(2,2) /jcob
      jcbi(1,2)=  - jcb(1,2) /jcob
      jcbi(2,1)=  - jcb(2,1) /jcob
      jcbi(2,2)=    jcb(1,1) /jcob

      do k=1,m
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
      end do

      do i=1,m
         i1=2*i-1
         i2=2*i
         Bmat(1,i1)=dNdx(i) ; Bmat(1,i2)=0.d0
         Bmat(2,i1)=0.d0    ; Bmat(2,i2)=dNdy(i)
         Bmat(3,i1)=dNdy(i) ; Bmat(3,i2)=dNdx(i)
      end do

      Ael=Ael + matmul(transpose(Bmat),matmul(penalty*Kmat,Bmat))*wq*jcob

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
   end do

   !==============================================!
   !=====[retrieve pressure]======================!
   !==============================================!

   do iel=1,nel

      rq=0.d0
      sq=0.d0
      
      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      jcbi(1,1)=    jcb(2,2) /jcob
      jcbi(1,2)=  - jcb(1,2) /jcob
      jcbi(2,1)=  - jcb(2,1) /jcob
      jcbi(2,2)=    jcb(1,1) /jcob

      do k=1,m
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
      end do

      xq=0.d0
      yq=0.d0
      exxq=0.d0
      eyyq=0.d0
      do k=1,m
         xq=xq+N(k)*x(icon(k,iel))
         yq=yq+N(k)*y(icon(k,iel))
         exxq=exxq+ dNdx(k)*u(icon(k,iel))
         eyyq=eyyq+ dNdy(k)*v(icon(k,iel))
      end do

      press(iel)=-penalty*(exxq+eyyq)

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

   if (iproc==0) call output_for_paraview2D (np,nel,x,y,u,v,press,icon)

   deallocate(x,y)
   deallocate(u,v)
   deallocate(icon)
   deallocate(bc_fix)
   deallocate(bc_val)
   deallocate(press)
   deallocate(idV%A_ELT)   
   deallocate(idV%RHS)    
   deallocate(idV%ELTPTR) 
   deallocate(idV%ELTVAR)     

end do ! kkk

if (iproc==0) close (1020)

end program

!==============================================!
!==============================================!
!==============================================!

function uth (x,y)
real(8) uth,x,y
uth = x**2 * (1.d0-x)**2 * (2.d0*y - 6.d0*y**2 + 4*y**3)
end function

function vth (x,y)
real(8) vth,x,y
vth = -y**2 * (1.d0-y)**2 * (2.d0*x - 6.d0*x**2 + 4*x**3)
end function

function pth (x,y)
real(8) pth,x,y
pth = x*(1.d0-x)
end function

function b1 (x,y)
real(8) b1,x,y
b1 = ( (12.d0-24.d0*y)*x**4 + (-24.d0+48.d0*y)*x**3 + (-48.d0*y+72.d0*y**2-48.d0*y**3+12.d0)*x**2 &
   + (-2.d0+24.d0*y-72.d0*y**2+48.d0*y**3)*x + 1.d0-4.d0*y+12.d0*y**2-8.d0*y**3 )
end function

function b2 (x,y)
real(8) b2,x,y
b2= ( (8.d0-48.d0*y+48.d0*y**2)*x**3 + (-12.d0+72.d0*y-72*y**2)*x**2 + &
    (4.d0-24.d0*y+48.d0*y**2-48.d0*y**3+24.d0*y**4)*x - 12.d0*y**2 + 24.d0*y**3 -12.d0*y**4)
end function



