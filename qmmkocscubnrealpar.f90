program MMK
    implicit none
    include "omp_lib.h"
    
    integer, parameter :: num_stream = 1
    integer, parameter :: N = 21, nt = 1 ! N need odd
    integer, parameter :: Nrealiz = 5
    complex(8), dimension(N, N) :: Hamiltonian, Ed ! , Hh, Hc, Hi
    complex(8), dimension(2*N, 2*N) :: c, Hcu, phprod, enprod, bb, AB, Hamilt , HttR
    complex(8), dimension(3,2*N) ::  Htt, Ht, Htt1, Htt2, Htt3, Htt4, Htt5
    complex(8), dimension(2*N) :: k1, k2, k3, k4, fosc, foc, focn, b, focn_m
    real(8), dimension(N) :: a1, a2, a
    real(8), dimension(2) :: fcub
    real, dimension(2) :: EigenValues
    complex, dimension(2,2) :: ff
    complex(8), dimension(2,2) :: Sx,  Sz , II, Splus, Sminus, S1, Sed
    integer :: i,j,k,p, tn, l, u, z, num
    real(8) :: s, ran, sumvec, ts, Pe, t, temp1, temp2, Posc, P1, P2
    integer, parameter :: nddt=30
    real(8), parameter :: pi = 3.14159265358979323846, dt = 0.1, ddt = dt/real(nddt)   
    real(8), parameter :: wq = 6. !qubit frequency
    real(8), parameter :: w = 6.001 !frequency of the field exciting the qubit
    real(8), parameter :: Gf_0 = 0. ! 0. !0.001 !0001 !0001 !phase relaxation
    real(8), parameter :: Ge_0 = 0. !0.00135 !00165 !165 ! 0.00075 !0.00065 !0.00165 !02 !02 !0.00065 !01 !energy relaxation
    real(8), parameter :: Go_0 = 0.01 !  01 !01 !0.0005 !1 !001 !001 !01 !01 ! 01 !0.01 !01 !02 !damping of the oscillator
    real(8) :: Ge, Go, Gf
    real(8), parameter :: wfo=1. !frequency of the field acting on the oscillator
    real(8), parameter :: f0= 0.1 !2 ! 2 !0.3 !amplitude of the oscillator field
    real(8), parameter :: Ams=0.02 !02 !02 !02 !qubit field amplitude
    real(8), parameter :: Mu  = 0. !0001 !0001 !0.0001 !001 !001 !non-linearity
    real(8), parameter :: lambda= 0.02 !02 !0.02 !02 !02 !02 !0.02 !Interaction
    real(8), parameter :: wj= 1.02 !1. !oscillator frequency
    real(8), parameter :: dwo=0. 
    real(8), parameter :: Wrab = sqrt(real((wq-w)*(wq-w)+Ams*Ams))
    real(8), parameter :: period = 2*pi/Wrab !Wrab
    integer, parameter :: tperiod = 1000, Nstep = int(2 * pi /dt), tmax = int(Nstep)*tperiod
    integer :: NumTime, countOsc
    real(8), dimension(tmax+1) :: qp1, nr, nr1, nr2, qp2, PeArr, norma !, nr1, nr2, nr, PeArr
    real(8), dimension(Nrealiz) :: nrnum
    real(8) :: t1, t2, Am, f, G
    logical :: flag1,  flag2, flag3, flag4, flag5, flag6, flag7
    complex(8), dimension(2*N,2*N) ::  Htmatrix
    

    open(unit=1,status='unknown',file='qp1.txt')
    open(unit=2,status='unknown',file='qp2.txt')
    open(unit=3,status='unknown',file='nr1.txt')
    open(unit=4,status='unknown',file='nr2.txt')
    open(unit=5,status='unknown',file='nr.txt')
    open(unit=7,status='unknown',file='nrnum.txt')
    
    
    
    
    
    flag1 = .false.
    flag2 = .false.
    flag3 = .false.
    flag4 = .false.
    flag5 = .false.
    flag6 = .false.

    Sx(1,1) = (0.,0.)
    Sx(1,2) = (1.,0.)
    Sx(2,1) = (1.,0.)
    Sx(2,2) = (0.,0.)

    
    Sz(1,1) = (1.,0.)
    Sz(1,2) = (0.,0.)
    Sz(2,1) = (0.,0.)
    Sz(2,2) = (-1.,0.)

    
    II(1,1) = (1.,0.)
    II(1,2) = (0.,0.)
    II(2,1) = (0.,0.)
    II(2,2) = (1.,0.)

    
    Splus(1,1) = (0.,0.)
    Splus(1,2) = (1.,0.)
    Splus(2,1) = (0.,0.)
    Splus(2,2) = (0.,0.)

    
    Sminus(1,1) = (0.,0.)
    Sminus(1,2) = (0.,0.)
    Sminus(2,1) = (1.,0.)
    Sminus(2,2) = (0.,0.)
    
    S1(1,1) = (1.,0.)
    S1(1,2) = (0.,0.)
    S1(2,1) = (0.,0.)
    S1(2,2) = (0.,0.)
    
    Sed(1,1) = (0.,0.)
    Sed(1,2) = (11.,0.)
    Sed(2,1) = (-11.,0.)
    Sed(2,2) = (1.,0.)
    
    
    
    
    do i = 1, N
   
        do j = 1, N
            Ed(i,j) = (0., 0.)
        enddo
        Ed(i,i) = (1., 0.)
    enddo
    

    do j = 1, tmax
        nr(j) = 0.
        nr1(j) = 0.
        nr2(j) = 0.
        qp1(j) = 0.
        qp2(j) = 0.
    enddo
    1. + Hc, 4)
    
    t1 = omp_get_wtime()
    call QuantumMethodMonteCarlo()
    t2 = omp_get_wtime()
    print *, "Time: ", t2 - t1, " sec    ", (t2 - t1)/60. ,"min"
    

    pause
    contains

    subroutine QuantumMethodMonteCarlo()
    integer :: p, i, j,z,l, k

            do p = 1, 2
                fcub(p)=(0.,0.) 
            end do
            do p = 1, 2*N
                fosc(p)=(0.,0.)
                foc(p)=(0,0.)
                focn(p)=(0,0.)
            end do
            fosc(nt)=(1.,0.)
            !qubit wave functions
            call EigenVecVal(H0(), EigenValues, ff, 2) 

            do p = 1, 2
                fcub(p)=ff(p, 1) 
            end do

            k=1
            do i = 1, N
                do j = 1, 2
                    foc(k) = fcub(j)*fosc(i)
                    k=k+1
                enddo
            enddo
            

            sumvec =0.
             do z = 1, 2*N
                sumvec = sumvec + abs(conjg(foc(z))*foc(z))
             end do 
             foc = foc/sqrt(sumvec)
            

    u=0
    num = 0
    countOsc=0
    !$omp parallel reduction(+:nr, qp1, nr1, nr2) private(Am, f, Ge, Gf, Go, Ht, Htt, flag1, flag2, flag3, flag4, flag5, temp1, Posc, P1, P2, i, j, k, l, z, t, ts, k1, k2, k3, k4, G, Pe, focn) shared(ran ,s, u) num_threads(num_stream)
    !$omp do
        do k = 1, Nrealiz
            flag1 = .false.
            flag2 = .false.
            flag3 = .false.
            flag4 = .false.
            flag5 = .false.
            
            focn=foc
            call random_seed()
            do j = 0, tmax	 
                t=j*dt 
                   
                 if (t<=1*period) then 
                    f = 0.
                    Am = Ams
			        Gf = 0.
			        Ge = 0.
			        Go = 0.
                    if (.NOT. flag1) then 
                        call formedStartHamiltInter(Htt)
                        flag1 = .true.
                        Ht = Htt
                    endif
                else
                    Go = Go_0
                    Am = 0.
                    if (.NOT. flag2) then 
                        call formedStartHamiltInter(Htt)   
                        flag2 = .true.
                        Ht = Htt
                    endif
                endif
                
                if (t>2*period) then 
                  f = f0
                  if (.NOT. flag4) then 
                    call formedStartHamiltInter(Htt)                     
                        flag4= .true.
                        Ht = Htt  
                  endif
                endif
                
                if (j>15000) then
                    Gf = Gf_0
                    Ge = Ge_0
                    if (.NOT. flag5) then 
                        call formedStartHamiltInter(Htt)
                            flag5= .true.
                            Ht = Htt
                     endif
                endif

                temp1 = 0.
                temp2 = 0.
                Posc = 0.
                P1 = 0.
                P2 = 0.
                do i = 1, N
                    temp1 = temp1 + abs(conjg(focn(2*i-1))*focn(2*i-1))
                    temp2 = temp2 + abs(conjg(focn(2*i))*focn(2*i)) 
                    Posc = Posc + (i-1)*(abs(conjg(focn(2*i-1))*focn(2*i-1)) + abs(conjg(focn(2*i))*focn(2*i)))
                    P1 = P1 + (i-1)*abs(conjg(focn(2*i-1))*focn(2*i-1))
                    P2 = P2 + (i-1)*abs(conjg(focn(2*i))*focn(2*i))
                enddo
                
                nr1(j+1) = nr1(j+1) + P1
                nr2(j+1) = nr2(j+1) + P2
                nr(j+1) = nr(j+1) + Posc
                qp1(j+1) = qp1(j+1) + temp1  
                qp2(j+1) = qp2(j+1) + temp2  
                
                
               Pe = 0.
               do l = 0, N-1
                    Pe = Pe + abs(conjg(focn(2*l+1))*focn(2*l+1))
               end do
               

                G = Gf + Ge*Pe+Go*Posc
                
                call RANDOM_NUMBER(ran)
                
                
                if (ran < dt*G) then
 
                  focn = focn/sqrt(sumvec)
                    call RANDOM_NUMBER(s)
                    if (s<Gf/G) then 
                        do l = 0, N-1
                            focn(2*l+2) = -focn(2*l+2)
                        enddo
   
                    else
                        if (s <(Gf/G + Ge*Pe/G)) then

                            do l = 1, N
                               focn(2*l)= focn(2*l-1)
                               focn(2*l-1)=(0.,0.)
                            enddo

                        else 

                            do l = 1, N-1
                                focn(2*l-1) = focn(2*l+1)*sqrt(real(l))
                                focn(2*l) = focn(2*l+2)*sqrt(real(l))
                            enddo
                            focn(2*N-1)= (0., 0.)
                            focn(2*N) = (0., 0.)
                            countOsc = countOsc + 1
                        endif
                    endif         
                else
                
                    ts = j*dt
                    do i=1, nddt
                        Ht(1,:) = Htt(1,:)*f*cos(wfo*ts)
                        Ht(2,:) = Htt(2,:)*Am*cos(w*ts)
                       k1 = ddt*RK(focn,ts,Ht)
                        Ht(1,:) = Htt(1,:)*f*cos(wfo*(ts+0.5*ddt))
                        Ht(2,:) = Htt(2,:)*Am*cos(w*(ts+0.5*ddt))
                       k2 = ddt*RK(focn+(0.5,0.)*k1,ts+0.5*ddt, Ht)
                       k3 = ddt*RK(focn+(0.5,0.)*k2,ts+0.5*ddt, Ht)
                        Ht(1,:) = Htt(1,:)*f*cos(wfo*(ts+ddt))
                        Ht(2,:) = Htt(2,:)*Am*cos(w*(ts+ddt))
                       k4 = ddt*RK(focn+k3,ts+ddt, Ht)
                       focn = focn + (1./6.)*(k1 + (2.,0.)*k2 + (2.,0.)*k3 + k4)
                       ts = ts + ddt  
                   end do
            endif

                  sumvec = 0.
                  do z = 1, 2*N
                      sumvec = sumvec + abs(conjg(focn(z))*focn(z))
                  end do 
                  focn = focn/sqrt(sumvec)
                  
            end do
            nrnum(k)=nr(tmax+1)
            u=u+1
            write (*,*) "Nrealiz", u
        end do
    !$OMP END DO    
    !$omp end parallel

        do i = 1, tmax
            t=(i-1)*dt
            WRITE(5,220) t/(2*pi), nr(i)/real(Nrealiz)
            WRITE(3,220) t/(2*pi), nr1(i)/real(Nrealiz)
            WRITE(4,220) t/(2*pi), nr2(i)/real(Nrealiz)
            WRITE(1,220) t/(2*pi), qp1(i)/real(Nrealiz)
            WRITE(2,220) t/(2*pi), qp2(i)/real(Nrealiz)
            if (i<=Nrealiz) then
                WRITE(7,221) nrnum(i)
            endif  
        end do
        
        
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        220 Format (F18.8, F18.8)
        221 Format (F18.8, F18.8)
        222 Format (F18.8)
        230 Format ('Progress: ', I3)
    end subroutine QuantumMethodMonteCarlo
    

    
    function RK(fs,t, H)
        complex(8), dimension(2*N) :: RK
        real(8), intent(in) :: t
        complex(8), dimension(2*N), intent(in) :: fs
        complex(8), parameter :: alpha = (0., -1.)
        complex(8), parameter :: beta = (0., 0.)
        complex(8), dimension(3,2*N), intent(in) :: H
        RK(1) = H(1,3)*fs(3)+H(2,2)*fs(2)+H(3,1)*fs(1)
        RK(2) = H(1,4)*fs(4)+H(2,3)*fs(3)+H(3,2)*fs(2)+H(2,2)*fs(1)
        do j = 1,2*N-4
            RK(j+2)=H(1,j+2)*fs(j)+H(2,j+2)*fs(j+1)+H(3,j+2)*fs(j+2)+H(2,j+3)*fs(j+3)+H(1,j+4)*fs(j+4)
        enddo

        RK(2*N-1) = H(2,2*N)*fs(2*N)+H(3,2*N-1)*fs(2*N-1)+H(2,2*N-1)*fs(2*N-2)+H(1,2*N-1)*fs(2*N-3)
        RK(2*N) = H(3,2*N)*fs(2*N)+H(2,2*N)*fs(2*N-1)+H(1,2*N)*fs(2*N-2)
        RK= alpha*RK
    end function RK
    
    
    function H0() 
        complex, dimension(2,2) :: H0
        H0(1,1) = (1.,0.)
        H0(1,2) = (0.,0.)
        H0(2,1) = (0.,0.)
        H0(2,2) = (-1.,0.)
        H0=(wq/2.)*H0
    end function H0
    
    subroutine formedStartHamiltInter(Htt)
        complex(8), dimension(3,2*N), intent(inout) :: Htt
        complex(8), dimension(2*N,2*N) ::  Hinter
        complex(8), dimension(2*N,2*N) ::  Hqubit
        complex(8), dimension(2*N,2*N) ::  Hoscill
        integer i, j, m
        complex(8), dimension(N,N) :: Hh, Hc, Hi
        do i = 1, N
            do j = 1, N
                Hc(i,j) = (0.,0.)
                Hh(i,j) = (0.,0.)
                Hi(i,j) = (0.,0.)
            end do
        end do

        do i = 1, 2*N
            do j = 1, 2*N
                Hoscill(i,j) = (0.,0.)
                Hinter(i,j) = (0.,0.)
                Hqubit(i,j) = (0.,0.)
                Htmatrix(i,j) = (0.,0.)
            end do
        end do
          
        do i = 0, N-1
            Hc(i+1,i+1) = i*wj - Mu*i*i - (0.,1.)*i*Go*0.5
        end do
        do i = 1, N-1
            Hh(i,i+1) = sqrt(real(i))
            Hh(i+1,i) = sqrt(real(i))
        enddo
        do i = 0, N-1
            Hi(i+1,i+1) = lambda*i*wj-Mu*i*i*lambda*0.25
        enddo
        do i = 1, 3
            do j = 1, 2*N
                Htt(i,j)=(0.,0.)
                Ht(i,j)=(0.,0.)
            end do
        end do

        
        call KronProd(Hi, Sz , Hinter)
        call KronProd(Ed, wq * Sz*0.5 + Sx*0.5 - (0.,1.)*Gf*0.5*II - (0.,1.)*Ge*0.5*S1 , Hqubit)
        call KronProd(Hh + Hc, II , Hoscill)
        Htmatrix = Hqubit + Hoscill + Hinter
        do j = 1, 2*N
            m=3-j
            do i = max(1, j - 2), j
                Htt(m + i, j) = Htmatrix(i, j)
            enddo
        enddo
    end subroutine formedStartHamiltInter
    
    subroutine KronProd(A,B,C)
       IMPLICIT NONE
       complex(8), dimension (:,:), intent(in)  :: A, B
       complex(8), dimension (size(A,1)*size(B,1),size(A,2)*size(B,2)), intent(inout) :: C
       integer :: i = 0, j = 0, k = 0, l = 0
       integer :: m = 0, nn = 0, p = 0, q = 0
       C=(0.,0.)
        do i = 1,size(A,1)
            do j = 1,size(A,2)
                nn=(i-1)*size(B,1) + 1
                m=nn+size(B,1) - 1
                p=(j-1)*size(B,2) + 1
                q=p+size(B,2) - 1
                C(nn:m,p:q) = A(i,j)*B
            enddo
        enddo 
    C=TRANSPOSE(C)
    end subroutine KronProd

    
    subroutine EigenVecVal(Matrix, EigenValues, EigenVectors, M)
        integer, intent(in) :: M
        complex, dimension(M,M), intent(inout) :: EigenVectors
        complex, dimension(M,M), intent(in) :: Matrix
        real, dimension(M), intent(inout) :: EigenValues
        integer, parameter :: LWMAX = 1000
        real, dimension(3*M-2) :: RWORK
        complex WORK(LWMAX)
        integer INFO, LWORK
        !documentation link:
        !software.intel.com/
        !cheev
        INFO=0
        LWORK = -1
        call cheev( 'V', 'U', M, Matrix, M, EigenValues, WORK, LWORK, RWORK, INFO)
        LWORK = MIN(LWMAX, INT(WORK(1)))
        call cheev( 'V', 'U', M, Matrix, M, EigenValues, WORK, LWORK, RWORK, INFO)
        EigenVectors = Matrix
        if( INFO.GT.0 ) then
            write(*,*)'The algorithm failed to compute eigenvalues.'
            stop
        end if
    end subroutine EigenVecVal
     
    subroutine PRINT_Eigenvalues(EigenValues, N)
        integer, intent(in) :: N
        real, dimension(1, *), intent(in) :: EigenValues
        integer :: i, j
        write(*,*) 'Eigenvalues'
        do i = 1, 1
            write(*,2) (EigenValues(i,j), j = 1, N)
        end do
        2 FORMAT( 11(:,1X,F15.6) )
    end subroutine PRINT_Eigenvalues
    
    end
