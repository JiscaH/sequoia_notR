! Calculate the expected distributions of:
! - differences between duplicate samples of same individual (dif_S), parent-offspring (dif_PO) and unrelated individuals (dif_U)
! - opposing homozygous loci between parent-offspring (OH_PO), full siblings (OH_FS), and unrelated (OH_U)
! - mendelian errors between parent-parent-offspring trio (ME_PPO)
!
! Jisca Huisman, 2022
!
! input: (text file)
! line 1: number of SNPs  (assuming missingness = 0!)
! line 2: mean MAF (or median)
! line 3-5: genotyping error matrix (3x3)
!
! output:
! line 1: header
! line 2: per-locus probability for each of these seven distributions
! line 3-(nq+2): quantiles from these seven distributions
! 
! Note: the quantiles are defined in 'program CalcMaxMismatch' below.
! when changing 'quants', adjust 'nq' (= length of quants) accordingly. 


!===============================================================================
!===============================================================================

module global
implicit none

  double precision :: ErrM(3,3), AHWE(3), OHWE(1,3), AKAP(3,3), AKA2P(3,3,3), &
    pOO_U(3,3), pOO_PO(3,3), pOO_FS(3,3)

contains
  
  !-----------------------------------------------
  subroutine init_inherit(MAF)
    implicit none
    
    double precision, intent(IN) :: MAF
    
    ! actual genotype frequencies
    AHWE(1) = (1-MAF)**2
    AHWE(2) = 2*MAF*(1-MAF)
    AHWE(3) = MAF**2
    
    ! observed genotype frequencies
    OHWE(1,:) = MATMUL(AHWE, ErrM)
    
    ! probability offspring genotype conditional on parent genotype, 
    ! other allele draw from population in HWE 
    AKAP(1, :) = (/ 1-MAF, (1-MAF)/2, 0.0D0 /)
    AKAP(2, :) = (/ MAF, 0.5D0, 1-MAF /)
    AKAP(3, :) = (/ 0D0, MAF/2, MAF /)
    
    ! inheritance probability conditional on both parents  (offspr - parent1 - parent2)
    AKA2P(1,1,:) = dble((/ 1.0, 0.5, 0.0 /))
    AKA2P(1,2,:) = dble((/ 0.5, 0.25, 0.0 /))
    AKA2P(1,3,:) = dble((/ 0.0, 0.0, 0.0 /))

    AKA2P(2,1,:) = dble((/ 0.0, 0.5, 1.0 /))
    AKA2P(2,2,:) = dble((/ 0.5, 0.5, 0.5 /))
    AKA2P(2,3,:) = dble((/ 1.0, 0.5, 0.0 /))

    AKA2P(3,1,:) = dble((/ 0.0, 0.0, 0.0 /))
    AKA2P(3,2,:) = dble((/ 0.0, 0.25, 0.5 /))
    AKA2P(3,3,:) = dble((/ 0.0, 0.5, 1.0 /))

  end subroutine init_inherit
  
  
  !-----------------------------------------------
  
  subroutine init_joined_distr()   ! joined distrbutions of observed genotypes
    implicit none
    
    double precision :: pAA_PO(3,3), pAA_FS(3,3), tmpM(3,3)
    integer :: i, j, g, h
    
    ! 2 unrelated individuals
    pOO_U = MATMUL(TRANSPOSE(OHWE), OHWE)
    
    ! parent-offspring
    do i=1, 3
      pAA_PO(i, :) = AKAP(i,:) * AHWE
    enddo
    pOO_PO = MATMUL( MATMUL( TRANSPOSE( ErrM ), pAA_PO ), ErrM )
    
    ! full siblings
    do g=1,3  ! sibling 1
      do h=1,3 ! sibling 2
        tmpM = 0D0
        do i=1,3  ! dam
          do j=1,3  ! sire 
            tmpM(i,j) = AKA2P(g,i,j) * AKA2P(h,i,j) * AHWE(i) * AHWE(j) 
          enddo
        enddo
        pAA_FS(g,h) = sum(tmpM)   ! sum over all possible parent genotypes
      enddo
    enddo
    pOO_FS = MATMUL( MATMUL( TRANSPOSE( ErrM ), pAA_FS ), ErrM )
  
  end subroutine init_joined_distr
  

  !-----------------------------------------------
  double precision function P_dif_S()
    implicit none

    double precision ::  tmpV(3)

    ! probability that genotype differs for duplicate sample = 1 - identical 
    ! = diagonal of ErrM weighed by genotype frequencies
    tmpV = MATMUL( ErrM**2, AHWE)
    P_dif_S = 1 - SUM(tmpV) 

  end function P_dif_S
  
  
  !-----------------------------------------------
  double precision function P_ME_PPO()
    implicit none

    double precision :: pAAA(3,3,3), pOOA(3,3,3), pOOO(3,3,3)
    integer :: i, j
    logical :: IME(3,3,3)
    
    ! indicator array: is mendelian error (1) or not (0)
    IME = .FALSE.
    IME(1,3,:) = .TRUE.
    IME(3,1,:) = .TRUE.
    IME(1,:,3) = .TRUE.
    IME(3,:,1) = .TRUE.
    IME(2,1,1) = .TRUE.
    IME(2,3,3) = .TRUE.
    
    do i=1,3
      do j=1,3
        pAAA(:,i,j) = AKA2P(:,i,j) * AHWE(i) * AHWE(j)
      enddo
    enddo
    
    do j=1,3
      pOOA(:,:,j) = MATMUL( MATMUL( TRANSPOSE( ErrM ), pAAA(:,:,j) ), ErrM )
    enddo
    
    do i=1,3
      pOOO(:,i,:) = MATMUL( pOOA(:,i,:), ErrM )
    enddo
    
    P_ME_PPO = SUM(pOOO, MASK = IME)

  end function P_ME_PPO

end module global


!===============================================================================

module binom
implicit none

integer, parameter :: max_fact_vals = 160   ! rounded to infinity for >170
double precision :: fact_vals(0:max_fact_vals)

contains  
  !-----------------------------------------------
  
  double precision function lfact(x)
    implicit none
    
    integer, intent(IN) :: x
    double precision :: z, pi
    
    z = 1.0D0 + x
    pi = 4 * DATAN(1.D0)  ! double precision arc-tan
    
    ! https://www.johndcook.com/blog/2010/08/16/how-to-compute-log-factorial/
    if (x <= max_fact_vals) then
      lfact = log( fact_vals(x) )
    else
      lfact = (z - 0.5) * log(z) - z + 0.5 * log(2*pi) + 1/(12*z)
    endif
  
  end function lfact
  
  !-----------------------------------------------
  
  double precision function ldbinom(x, n, p)
    implicit none

    integer, intent(IN) :: x, n
    double precision, intent(IN) :: p
    
    ldbinom = lfact(n) - lfact(x) - lfact(n-x) + x*log(p) + (n-x)*log(1-p)

  end function ldbinom

  !-----------------------------------------------

  integer function qbinom(q, n, p)
    implicit none

    integer, intent(IN) :: n
    double precision, intent(IN) :: q, p
    double precision:: probs(0:n), cumprobs(0:n)
    integer :: i

    if (abs(q - 0.0D0) < 1e-16) then
      qbinom = 0
      
    else if (abs(q - 1.0D0) < 1e-16) then
      qbinom = n
    
    else
      probs = 0D0
      cumprobs = 0D0
      do i = 0, n
        probs(i) = ldbinom(i, n, p)
        cumprobs(i) = SUM( exp( probs(0:i) ) )
      enddo
      
      qbinom = MINLOC(cumprobs, MASK = cumprobs >= q, DIM = 1) -1   ! starts at 0
    endif

  end function qbinom
  
  !-----------------------------------------------
  
  subroutine calc_fact_vals
    implicit none
    
    integer :: i
      
    fact_vals(0) = 1.0D0
    do i = 1, max_fact_vals
      fact_vals(i) = i * fact_vals(i-1)
    enddo

  end subroutine calc_fact_vals

  
end module binom



!===============================================================================
!===============================================================================

program CalcMaxMismatch
  use global
  use binom
  implicit none

  integer :: nSnps, i, x
  integer, parameter :: nq = 11  ! number of quantiles
  double precision :: MAF,  quants(nq), pOUT(7)
  integer :: nOUT(nq, 7)
  
  quants = DBLE( (/ 0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 0.9999, 0.99999, 1.0 /) )


  call calc_fact_vals

  call ReadInput(nSnps, MAF, ErrM)
  
  call init_inherit(MAF)

  call init_joined_distr()

  pOUT = 0D0
  pOUT(1) = P_dif_S()
  pOUT(2) = 1 - pOO_PO(1,1) - pOO_PO(2,2) - pOO_PO(3,3)
  pOUT(3) = 1 - pOO_U(1,1) - pOO_U(2,2) - pOO_U(3,3)
  pOUT(4) = pOO_PO(1,3) + pOO_PO(3,1)
  pOUT(5) = pOO_FS(1,3) + pOO_FS(3,1)
  pOUT(6) = pOO_U(1,3) + pOO_U(3,1)
  pOUT(7) = P_ME_PPO()

  do x= 1, 7
    do i=1, nq
      nOUT(i,x) = qbinom(quants(i), nSnps, pOUT(x))
    enddo
  enddo


  ! write output
  open (unit = 201, file = 'Mismatch_distr.txt', status = 'unknown')
    write (201, '(a8, 2X, 7a8)')  'quantile', 'dif_S', 'dif_PO', 'dif_U', 'OH_PO', 'OH_FS', 'OH_U', 'ME_PPO'
    write(201, '(4X, "prob", 2X, 7f8.4)') pOUT(:)
    do i = 1, nq
      write(201, '(f8.6, 2X, 7i8)')  quants(i), nOUT(i,:)
    enddo
  close(201)
  
  print *, 'Done.'

end program CalcMaxMismatch


!===============================================================================
!===============================================================================

subroutine ReadInput(nSnps, MAF, ErrM)
  implicit none

  integer, intent(OUT) :: nSnps
  double precision, intent(OUT) :: MAF, ErrM(3,3)
  character(len = 200) :: tag
  integer :: i

  nSnps = 0
  MAF = 0D0
  ErrM = 0D0

  open (unit=101, file="MismatchSpecs.txt", status="old")
  read(101, *) tag, nSnps
  read(101, *) tag, MAF
  read(101, *) tag, ErrM(1,1:3)
  read(101, *) tag, ErrM(2,1:3)
  read(101, *) tag, ErrM(3,1:3)
  close(101)
  
  do i=1,3
    if (abs( sum(ErrM(i,:)) - 1.0D0 ) > 0.0001)  stop('rows in ErrM must sum to 1')
  enddo
  

end subroutine ReadInput

!===============================================================================