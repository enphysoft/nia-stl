module math
   implicit none
   public :: rkind, ikind, ckind !, success, i,j,k, ip, jp, kp,ihalf, xrand 
   integer   ,parameter :: rdigits=8
   integer   ,parameter :: idigits=8
   integer   ,parameter :: idigits4=4
   integer   ,parameter :: rkind=8, ikind=4, ikind4=4, ckind=1

   real(kind=rkind)  ,parameter :: pi_val=3.141592653589793238462643383279
   real(kind=rkind)  ,parameter :: ZERO=0.0D0 
   real(kind=rkind)  ,parameter :: ONE=1.0D0   ,TWO=2.0D0   ,THREE=3.0D0
   real(kind=rkind)  ,parameter :: FOUR=4.0D0  ,FIVE=5.0D0  ,SIX=6.0D0
   real(kind=rkind)  ,parameter :: SEVEN=7.0D0 ,EIGHT=8.0D0 ,NINE=9.0D0, TEN=10.0D0

   real(kind=rkind)  :: nmvecX(3), nmvecY(3), nmvecZ(3)  
   data nmvecX /1.0d0, 0.0d0, 0.0d0/
   data nmvecY /0.0d0, 1.0d0, 0.0d0/ 
   data nmvecZ /0.0d0, 0.0d0, 1.0d0/

   interface operator (.norm.)
     module procedure OP_NORM
   end interface
   
   interface operator (.dot.)
     module procedure OP_DOT
   end interface
   
   interface operator (.cross.)
     module procedure OP_CROSS
   end interface

contains

   integer function iself0 (i) result(j)
      implicit none
      integer, intent(inout) :: i
      j = i 
   end function iself0
      
   integer function imnus1 (i) result(j)
      implicit none
      integer, intent(inout) :: i
      j = i - 1
      i = j
   end function imnus1
   
   integer function iplus1 (i) result(j)
      implicit none
      integer, intent(inout) :: i
      j = i + 1
      i = j
   end function iplus1
      
   logical function cross(Avec,Bvec,Cvec) result(message)
   ! A x B = C
     real (rkind) :: Avec(3), Bvec(3), Cvec(3)
     Cvec(1)=Avec(2)*Bvec(3)-Bvec(2)*Avec(3)
     Cvec(2)=Avec(3)*Bvec(1)-Bvec(3)*Avec(1)
     Cvec(3)=Avec(1)*Bvec(2)-Bvec(1)*Avec(2)
     message=.true.
   end function cross
   
   function OP_NORM (A,B) result (C)
     implicit none
     real (rkind), dimension(:), intent(in)    :: A, B
     real (rkind) :: C
     integer :: na, nb, i
     na = SIZE(A)
     nb = SIZE(B)
     C  = OP_DOT (A,B)
     C  = SQRT(C)
   end function OP_NORM
   
   function OP_DOT (A,B) result (C)
   ! A dot B = C
     implicit none
     real (rkind), dimension(:), intent(in)    :: A, B
     real (rkind) :: C
     integer :: na, nb, i
     na = SIZE(A)
     nb = SIZE(B)
     if (na .eq. nb) then
        C = 0.0d0
        do i = 1, na
          C = C + A(i) * B(i)
        enddo
     else
     endif
   end function OP_DOT
   
   function OP_CROSS (A,B) result (C)
   ! A x B = C
     implicit none
     real (rkind), dimension(3), intent(in)    :: A, B
     real (rkind), dimension(3)     :: C
     C(1) = A(2)*B(3) - A(3)* B(2)
     C(2) = A(3)*B(1) - A(1)* B(3)
     C(3) = A(1)*B(2) - A(2)* B(1)
     return
   end function OP_CROSS
   
   function magntd_vector(A) result (magntd)
     implicit none
     real (rkind), dimension(3), intent(in)    :: A
     real (rkind) :: magntd
     magntd = dsqrt(A .dot. A)
   
   end function
   
   function uvec(A) result (Aunit)
   ! This is a masking fuction of unit_vector(A) to shorten the name.
     implicit none
     real (rkind), dimension(3), intent(in)    :: A
     real (rkind), dimension(3) :: Aunit
     Aunit = unit_vector(A)
   end function uvec
   
   
   function unit_vector(A) result (Aunit)
     implicit none
     real (rkind), dimension(3), intent(in)    :: A
     real (rkind), dimension(3) :: Aunit
     real (rkind) :: magntd
   !  magntd = dsqrt(A .dot. A)
     magntd = magntd_vector(A)  
     if(magntd.eq.0.0d0) then 
       Aunit = magntd
     else
       Aunit = A / magntd
     endif
     return
   end function unit_vector
   
   subroutine vector_renorm (vector, magntd, unit_vector)
     implicit none
     real (rkind), dimension(3), intent(in)    :: vector
     real (rkind), dimension(3), intent(out)    :: unit_vector
     real (rkind), intent(out)    :: magntd
     magntd = dsqrt(vector .dot. vector  )
     if(magntd.eq.0.0d0) then
       unit_vector = 0.0d0
     else
       unit_vector = vector / magntd
     endif
     return
   end subroutine vector_renorm
   
   subroutine ortho_transform (x1,y1, theta, x2,y2)
     implicit none
     real (rkind),  intent(in)    :: x1, y1, theta
     real (rkind),  intent(out)   :: x2, y2
     real (rkind),  dimension(2,2):: P
     logical                      :: message
     message =  polar_matrix (theta, P)
     x2 = P(1,1) * x1 + P(1,2)*y1
     y2 = P(2,1) * x1 + P(2,2)*y1
     return
   end subroutine ortho_transform
   
   
   logical function polar_matrix (theta, polarMat) result(message)
     implicit none
     real (rkind),  intent(in)    :: theta
     real (rkind),  dimension(2,2),  intent(out)   :: polarMat
   !  x' =   cos(t)*x + sin(t)*y
   !  y' =  -sin(t)*x + cos(t)*y
     polarMat(1,1) =   DCOS(theta)   ;   polarMat(1,2) =   DSIN(theta)
     polarMat(2,1) = - polarMat(1,2) ;   polarMat(2,2) = + polarMat(1,1)
     message = .true.
     return
   end function polar_matrix
   
end module math
 
! EOF    
