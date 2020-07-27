module mathstl
 ! useful web site for trangular math
   ! http://mathworld.wolfram.com/TriangleArea.html
   ! http://mathworld.wolfram.com/Line-PlaneIntersection.html
   use math
   use fufs  
   type vector 
      real(8) :: postn(3)
   end type vector

   type triVertex
      ! primary 
      real(8)   :: vertexMat(3,3)       ! each row for vertexA, B, and C. 
      real(8)   :: vertexA(3)
      real(8)   :: vertexB(3)
      real(8)   :: vertexC(3)
      real(8)   :: vertexM(3)     !to replace tmid 
      ! derivative
      integer   :: zoneidA(0:2)         ! for each of vertex A, B, and C 
      integer   :: zoneidB(0:2)         !      (ipx, ipy) => (irx, iry) 
      integer   :: zoneidC(0:2)         !      => address = Nrx * (iry-1) + irx 
      integer   :: zoneidM(0:2)
   end type triVertex
   
   type triVertexMat
      real(8)   :: vertexMat(3,3) ! each row for vertexA, B, and C. 
      integer   :: regionID(3)    ! ╚ following the structure of STL facet component
   end type triVertexMat
   
   type relTriVectors
      real(8) :: vecBwrtA(3)    ! with respect to the previous one
      real(8) :: vecCwrtB(3)
      real(8) :: vecAwrtC(3)
      real(8) :: vecAwrtM(3)    ! with respect to the mid centroid.
      real(8) :: vecBwrtM(3)
      real(8) :: vecCwrtM(3)
   end type relTriVectors

   type stl_list
      integer :: ix
      integer :: jy
   end type stl_list

   type stl_tri
      integer                :: tid          ! triangle ID, sequential number
      integer                :: cid          ! coordinate normal vector id: 1,2,3 = x,y,z, unless 0
      integer                :: lat_nnbor(3) ! lateral-sharing neighbors ID 1,2,3
      integer                :: vtx_nnbor(8) ! vortex-sharing neighbors
      character(1)           :: bid          ! L, R, F, B - boundary association character 
      real(8)                :: area         ! triangular area 
      real(8)                :: tnvec(3)     ! normal vector
      real(8)                :: tmid (3)     ! midpoint = centroid to be calculated 
      type(triVertex)        :: tvtx         ! vertex points: vertexA, B, and C.
      type(relTriVectors)    :: trvs         ! triangle relative vectors 
   end type stl_tri
   
   type stl_nnb
      integer :: fidPair(3)
      integer :: vidPair(30)
      integer :: vidPairMax
      integer :: matPair(3,3)
      real(8) :: mvtx(3,3)   ! mvtx(i,j) i_th vertex's k_th component
   end type stl_nnb
   
   type stl_tri_nnb
      type (stl_tri) :: fct
      type (stl_nnb) :: nnb      
   end type stl_tri_nnb
   
contains
   
   !
   ! LIST of functions and 
   ! 
   ! DOUBLE PRECISION function triArea (postnA,postnB,postnC) RESULT (Area)
   ! DOUBLE PRECISION function M44DET (A) RESULT (DET)
   ! DOUBLE PRECISION function M33DET (A) RESULT (DET)
   ! integer function check_uvec_dvec (uvec)  result(cid)
   ! DOUBLE PRECISION function stl_interpol_z (hx1)  result (zpostn)

   !****************************************************************
   !
   !                            functionS
   !
   !****************************************************************
   logical function chk_equal_vectors (vectA, vectB) result (chk_message)
      ! [2020-01-05-13-19-18-PM] 
      real (8), intent (in) :: vectA(3), vectB(3)
      chk_message = .false.
      if ( (vectA(1) == vectB(1))  .and. &
           (vectA(2) == vectB(2))  .and. &
           (vectA(3) == vectB(3))) then
         chk_message = .true.
      end if
   end function chk_equal_vectors



   integer function cal_num_of_file_lines (inputfile) result (numline)
      ! [2020-01-05-13-19-30-PM] 
      implicit none
      ! This function is to count the number of lines of inputfile.
      ! This works with only gfortran, and ifort give erroneous number.
      character(128), intent(in)   :: inputfile
      character(128)               :: inputline
      integer                      :: eastat, lin, numvalues
      lin = 10
      eastat = 0 
      numvalues = 0 
      
      OPEN(lin, file=inputfile , status='old')
      loop1: DO
         READ(lin,*,iostat=eastat) inputline
         ! write(*,*)  "eastat = ", eastat
         numvalues = numvalues + 1
         IF (eastat < 0) THEN
            WRITE(*,*)  'Number of lines in ', trim(inputfile),": ",numvalues-1
            numline   = numvalues - 1
            EXIT loop1
         ELSE IF (eastat > 0) THEN
            write(*,* ) "IOSTAT is ", eastat
            STOP 'IO-error'
         ENDIF
      END DO loop1
      CLOSE(lin)
   end function cal_num_of_file_lines

   !****************************************************************
   ! triArea - calculates the triangle area 
   !*****************************************************
   real*8 function triArea (postnA,postnB,postnC) RESULT (Area)
      implicit none
      real (8)  :: postnA(3),postnB(3),postnC(3)
      real (8)  :: vecB_A(3), vecC_A(3), crossVec(3)
      vecB_A = postnB - postnA
      vecC_A = postnC - postnA
      crossVec = vecB_A  .cross. vecC_A
      area     = sqrt(abs (crossVec .dot. crossVec)) / 2.0d0
      return 
   END function triArea
   !****************************************************************
   ! triArea2D - to calculate an area of 2D triangle in x-y plane
   !****************************************************************
   real*8 function triArea2D (postn2d_A,postn2d_B,postn2d_C) RESULT (Area2d)
      IMPLICIT NONE
      real*8 , intent(in)       :: postn2d_A(2), postn2d_B(2), postn2d_C(2)
      real*8                    :: triMat(3,3), detTriMat
      TriMat (:,:)   = 1.0d0
      TriMat (1,1:2) = postn2d_A
      TriMat (2,1:2) = postn2d_B
      TriMat (3,1:2) = postn2d_C
      detTriMat         = M33DET (TriMat)
      Area2D            = 0.50d0* abs(detTriMat)
   end function triArea2D
   !*****************************************************
   !  M44DET  -  Compute the determinant of a 4x4 matrix.
   !*****************************************************
   real*8 function M44DET (A) RESULT (DET)
      IMPLICIT NONE
      real*8, DIMENSION(4,4), INTENT(IN)  :: A
      DET= A(1,1)*(&
           A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+&
           A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+&
           A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))-&
           A(1,2)*(&
           A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+&
           A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
           A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+&
           A(1,3)*(&
           A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+&
           A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
           A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-&
           A(1,4)*(&
           A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+&
           A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+&
           A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
      RETURN
   END function M44DET
   !******************************************************
   !  M33DET  -  Compute the determinant of a 3x3 matrix.
   !******************************************************
   real*8 function M33DET (A) RESULT (DET)
      IMPLICIT NONE
      real*8, DIMENSION(3,3), INTENT(IN) :: A
      DET =  A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  &
           - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)  &
           + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
      RETURN
   END function M33DET
   !****************************************************************
   !  check_uvec_dvec - check one of unit vector components is 1.0
   !****************************************************************
   integer function check_uvec_dvec (uvec)  result(cid)
      implicit none  
      real(8), intent(in)       :: uvec(3)
      integer   :: i 
      cid       = 0
      do i = 1, 3
         if (abs(uvec(i)) == 1.0D0)  then
            cid = i
         endif
      enddo
   END function check_uvec_dvec
   !****************************************************************
   !  stl_interpol_z - calculate the interpolated position
   !****************************************************************
   real*8 function stl_interpol_z (hx1)  result (zpostn)
      implicit none
      type(stl_tri), intent(in) :: hx1
      real(8) :: mat44a(4,4), mat44b(4,4), det44a, det44b
      mat44a(1,:)    = 1.0D0
      mat44a(2:4,1)  = hx1%tvtx%vertexA
      mat44a(2:4,2)  = hx1%tvtx%vertexB
      mat44a(2:4,3)  = hx1%tvtx%vertexC
      mat44a(2:4,4)  = hx1%tmid
      mat44a(4,4)    = 0.0d0
      mat44b         = mat44a
      mat44b(:,4)    = 0.0d0
      mat44b(4,4)    = 1.0d0
      det44a         = M44DET (mat44a)
      det44b         = M44DET (mat44b)
      zpostn  = - det44a/det44b
   END function stl_interpol_z
   !*********************************************************************
   !  check_point_inTri - to check a point is within a triangle.
   !*********************************************************************
   logical function check_point_inTri_Area  &
         (interpostn,hx1,nout,verbose) result (check_point_inTri_message)
      implicit none
      real(8)           , intent(in) :: interpostn(3)
      type(stl_tri)     , intent(in) :: hx1
      real(8)   :: vecIwrtA (3),vecIwrtB(3),vecIwrtC(3), val_loc_test(3)
      real(8)   :: normalVec(3), pArea(3), sumArea
      integer   :: i, itriSign(3),nout
      logical   :: verbose
      check_point_inTri_message = .true. 
      call calc_tri_area (interpostn,hx1%tvtx%vertexA,hx1%tvtx%vertexB,pArea(1))
      call calc_tri_area (interpostn,hx1%tvtx%vertexB,hx1%tvtx%vertexC,pArea(2))
      call calc_tri_area (interpostn,hx1%tvtx%vertexC,hx1%tvtx%vertexA,pArea(3))
      sumArea = pArea(1) + pArea(2) +  pArea(3) 
      sumArea = sumArea * (1.0D0)
      if (sumArea .gt. hx1%area ) check_point_inTri_message = .false.
      if (verbose .eqv. .true.) &
           write(nout,*) "hx1%area, sumArea", check_point_inTri_message, hx1%area, sumArea
   END function check_point_inTri_Area
   !*********************************************************************
   !  check_point_inTri - to check a point is within a triangle.
   !*********************************************************************
   logical function check_point_inTri  &
        (interpostn,hx1,verbose) result (check_point_inTri_message)
      implicit none
      real(8)           , intent(in) :: interpostn(3)
      type(stl_tri)     , intent(in) :: hx1
      real(8)   :: vecIwrtA (3),vecIwrtB(3),vecIwrtC(3), val_loc_test(3)
      real(8)   :: normalVec(3)
      integer   :: i, itriSign(3)
      logical   :: verbose
      vecIwrtA     = interpostn -  hx1%tvtx%vertexA
      vecIwrtB     = interpostn -  hx1%tvtx%vertexB
      vecIwrtC     = interpostn -  hx1%tvtx%vertexC
      normalVec    =  (hx1%trvs%vecBwrtA .cross. hx1%trvs%vecAwrtC)
      val_loc_test(1) = 1.0d0
      val_loc_test(1) = (vecIwrtA .cross. hx1%trvs%vecBwrtA) .dot. normalVec
      val_loc_test(2) = (vecIwrtB .cross. hx1%trvs%vecCwrtB) .dot. normalVec
      val_loc_test(3) = (vecIwrtC .cross. hx1%trvs%vecAwrtC) .dot. normalVec
      check_point_inTri_message = .true. 
      ! do i = 1, 3 
      !    if (val_loc_test(i) .lt. 0.0) then
      !       check_point_inTri_message = .false.
      !       exit
      !    endif
      ! enddo
      do i = 1, 3 
         itriSign(i)= sign(val_loc_test(i),abs(val_loc_test(i)))
      enddo 
      if ( itriSign(1) == itriSign(2) .and. itriSign(2) == itriSign(3)) then
         check_point_inTri_message = .true.
      else
         check_point_inTri_message = .false.
      endif
      if (verbose .eqv. .true.)  then
         write(*,"('check_point_inTri_message = ',L1,': ',4(2x,F18.9))") & 
              check_point_inTri_message, val_loc_test
      end if
   END function check_point_inTri
!****************************************************************
! 
! SUBROUTINES
! 
!****************************************************************
   !****************************************************************
   ! calc_tri_area - to calc. triangular area in 3D
   !****************************************************************
   subroutine CALC_TRI_AREA (postnA,postnB,postnC,triArea)
      implicit none
      real (8)  :: postnA(3),postnB(3),postnC(3)
      real (8)  :: vecB_A(3), vecC_A(3), crossVec(3),triArea
      vecB_A    = postnB - postnA
      vecC_A    = postnC - postnA
      crossVec  = vecB_A  .cross. vecC_A
      triArea   = abs (sqrt(crossVec .dot. crossVec)) / 2.0D0
   end subroutine CALC_TRI_AREA
   !****************************************************************
   ! calc_tri_area - to calc. triangular area in 3D
   !****************************************************************
   subroutine CALC_TRI_AREA2D (postn2d_A,postn2d_B,postn2d_C,triArea2d)
      implicit none
      real (8) :: postn2d_A(2), postn2d_B(2), postn2d_C(2)
      real (8) :: triMat(3,3), detTriMat
      real (8) :: triArea2d
      TriMat (:,:)   = 1.0d0
      TriMat (1,1:2) = postn2d_A
      TriMat (2,1:2) = postn2d_B
      TriMat (3,1:2) = postn2d_C
      detTriMat         = M33DET (TriMat)
      triArea2d         = 0.50d0* abs(detTriMat)
   end subroutine CALC_TRI_AREA2D
   !****************************************************************
   ! calc_tri_area - to calc. triangular area
   !****************************************************************
   subroutine CALC_PARALLELOGRAM_AREA_SQUARE (postnA,postnB,postnC,paraAreaSquare)
      implicit none
      real (8)  :: postnA(3),postnB(3),postnC(3)
      real (8)  :: vecB_A(3), vecC_A(3), crossVec(3),paraAreaSquare
      vecB_A    = postnB - postnA
      vecC_A    = postnC - postnA
      crossVec  = vecB_A  .cross. vecC_A
      paraAreaSquare  = crossVec .dot. crossVec
   end subroutine CALC_PARALLELOGRAM_AREA_SQUARE
   !****************************************************************
   !  stl_interpol_z - calculate the interpolated position, newest 
   !****************************************************************
   subroutine calc_stl_interpol_zpostn (itri,hx1,interpostn,zpostn,nout,verbose)
      ! This is to calculate the interpolated position of a triangle.
      ! 2018-07-19-08-16-42-AM 
      implicit none
      integer,          intent(in) :: itri 
      type(stl_tri),    intent(in) :: hx1
      logical                   :: verbose
      real(8), intent(inout)    :: interpostn(3)
      real(8)   :: mat44a(4,4), mat44b(4,4), det44a, det44b, zpostn
      integer   :: nout
      mat44a(1,:)       = 1.0D0
      mat44a(2:4,1)     = hx1%tvtx%vertexA
      mat44a(2:4,2)     = hx1%tvtx%vertexB
      mat44a(2:4,3)     = hx1%tvtx%vertexC
      mat44a(2:4,4)     = interpostn
      ! mat44a(4,4)    = 0.0d0
      mat44b            = mat44a
      mat44b(:,4)       = 0.0d0
      mat44b(4,4)       = 1.0d0
      det44a            = M44DET (mat44a)
      det44b            = M44DET (mat44b)
      if (det44b .ne. 0.0d0) then
         zpostn            = - det44a/det44b
         interpostn(3)     = interpostn(3) + zpostn
      else
         interpostn(3)     = interpostn(3) + hx1%tmid(3)
         write(900,*) "det44b is zero in subroutine CALC_STL_INTERPOL_ZPOSTN in itri = " , itri 
         write(900,*) "interpostn(3)     = interpostn(3) + hx1%tmid", interpostn(3)
         write(900,*) "DET mat44a(,)", det44a
         write(900,"('A=[',3(2X,F25.16,','),(2X,F25.16),';')") mat44a(1,:)
         write(900,"(      3(2X,F25.16,','),(2X,F25.16),';')") mat44a(2,:)
         write(900,"(      3(2X,F25.16,','),(2X,F25.16),';')") mat44a(3,:)
         write(900,"(      3(2X,F25.16,','),(2X,F25.16),']')") mat44a(4,:)
         write(900,*) "DET mat44b(,)", det44b
         write(900,"('B=[',3(2X,F25.16,','),(2X,F25.16),';')") mat44b(1,:)
         write(900,"(      3(2X,F25.16,','),(2X,F25.16),';')") mat44b(2,:)
         write(900,"(      3(2X,F25.16,','),(2X,F25.16),';')") mat44b(3,:)
         write(900,"(      3(2X,F25.16,','),(2X,F25.16),']')") mat44b(4,:)
         write(900,*) 
      endif
      if (verbose .eqv. .true.) &
           call print_info (itri,hx1,mat44a,mat44b,det44a,det44b,zpostn,nout,verbose)
      return
   end subroutine CALC_STL_INTERPOL_ZPOSTN
   !****************************************************************
   !  stl_interpol_z - calculate the interpolated position
   !****************************************************************
   subroutine CALC_STL_INTERPOL_Z (itri,hx1,zpostn,verbose)
      ! This is to calculate the interpolated position of a triangle.
      ! 2018-07-19-08-16-42-AM 
      implicit none
      integer,          intent(in) :: itri 
      type(stl_tri),    intent(in) :: hx1
      logical                   :: verbose
      real(8)   :: mat44a(4,4), mat44b(4,4), det44a, det44b, zpostn
      integer   :: nout=6
      mat44a(1,:)    = 1.0D0
      mat44a(2:4,1)  = hx1%tvtx%vertexA
      mat44a(2:4,2)  = hx1%tvtx%vertexB
      mat44a(2:4,3)  = hx1%tvtx%vertexC
      mat44a(2:4,4)  = hx1%tmid
      mat44a(4,4)    = 0.0d0
      mat44b         = mat44a
      mat44b(:,4)    = 0.0d0
      mat44b(4,4)    = 1.0d0
      det44a         = M44DET (mat44a)
      det44b         = M44DET (mat44b)
      zpostn  = - det44a/det44b
      if (verbose .eqv. .true.) & 
           call print_info (itri,hx1,mat44a,mat44b,det44a,det44b,zpostn,nout,verbose)
      return
   end subroutine CALC_STL_INTERPOL_Z
   !*********************************************************************
   !  analyze_ascii_stl - to analyze structural properties of a triangle
   !*********************************************************************
   subroutine ANALYZE_ALL_TRIANGLES (ntri,hx)
      implicit none 
      integer      , intent(in)         :: ntri
      type(stl_tri), intent(inout)      :: hx(ntri)
      integer                           :: itri
      do itri = 1, ntri
         call ANALYZE_A_TRIANGLE (hx(itri))
      enddo
   end subroutine ANALYZE_ALL_TRIANGLES
   !*********************************************************************
   !  analyze_ascii_stl - to analyze structural properties of a triangle
   !*********************************************************************
   subroutine ANALYZE_A_TRIANGLE (hx1) ! for cid, tmid, relpostn_ABC, area
      implicit none 
      type(stl_tri),            intent(inout)   :: hx1
      integer                   :: cid
      real (8)  :: crossVector(3), triArea, mat33(3,3)
      cid       = check_uvec_dvec (hx1%tnvec)
      hx1%cid   = cid 
      hx1%tvtx%vertexM  = (     hx1%tvtx%vertexA + & 
                                hx1%tvtx%vertexB + & 
                                hx1%tvtx%vertexC        ) / 3.0D0
      hx1%tmid          = hx1%tvtx%vertexM 
      hx1%trvs%vecBwrtA = hx1%tvtx%vertexB - hx1%tvtx%vertexA 
      hx1%trvs%vecCwrtB = hx1%tvtx%vertexC - hx1%tvtx%vertexB
      hx1%trvs%vecAwrtC = hx1%tvtx%vertexA - hx1%tvtx%vertexC
      mat33(1,:)        = hx1%tvtx%vertexA
      mat33(2,:)        = hx1%tvtx%vertexB
      mat33(3,:)        = hx1%tvtx%vertexC
      hx1%tvtx%vertexMat= mat33 ! Each row (of 3 element) is a vertex position. 
      call calc_tri_area (hx1%tvtx%vertexA, hx1%tvtx%vertexB,hx1%tvtx%vertexC,triArea)
      hx1%area  = triArea 
   end subroutine ANALYZE_A_TRIANGLE
   !*********************************************************************
   !  make_xy_grid - to make xy grid for interpolation  
   !*********************************************************************
   subroutine MAKE_XY_GRID  & 
        (Nx,Ny,Min_X,Min_Y,Max_X,Max_Y,dx,dy,marginX,marginY,Length,Width,postnX,postnY)
      use fufs
      implicit none
      integer , intent(in)      :: Nx,Ny
      real(8) , intent(in)      :: Min_X,Min_Y,Max_X,Max_Y,marginX,marginY,Length,Width,dx,dy
      real(8) , intent(inout)   :: postnX(Nx),postnY(Ny)
      integer :: i, j , nxygrid=10 
      real(8) :: Xpostn_Ref, Ypostn_Ref
      Xpostn_Ref = Min_X + marginX
      Ypostn_Ref = Min_Y + marginY
      Do i = 1, Nx
         postnX (i) = Xpostn_Ref + dx * dble(i-1) ; ! write(*,*)  postnX (i)      
      enddo
      Do j = 1, Ny
         postnY (j) = Ypostn_Ref + dy * dble(j-1) ; ! write(*,*)  postnY (j)  
      enddo
      ! 
   end subroutine MAKE_XY_GRID
   !*********************************************************************
   ! make_ray_vector - to generate a ray vector using three points, (x1, y1, z1)
   ! 2018-07-20-21-53-45-PM
   !*********************************************************************
   subroutine MAKE_RAY_VECTOR (ray_vector, x1,y1,z1)
      implicit none
      real(8), intent(in)       :: x1, y1, z1
      real(8), intent(out)      :: ray_vector(3)
      ray_vector(1) = x1
      ray_vector(2) = y1
      ray_vector(3) = z1
      return
   end subroutine MAKE_RAY_VECTOR
   !*********************************************************************
   !  check_ray_triangle_intersection - to check a point is within a triangle.
   !*********************************************************************
   subroutine check_point_triangle_intersection_2d (rayvec, hx1, nout, verbose, message)
      ! This routine calculates faster than "check_ray_triangle_intersection", 
      ! proviing an idential result. 
      ! 
      implicit none
      real(8)           , intent(in)    :: rayvec(3)
      type(stl_tri)     , intent(in)    :: hx1
      integer           , intent(in)    :: nout
      logical           , intent(in)    :: verbose
      logical           , intent(inout) :: message
      real(8)   :: normalVec(3), pArea(3), sumArea, tiny=1.E-3, nzmag
      real(8)   :: pseudoParaAreaSquare, tinyNum=1.0E-3, hxArea2D
      nzmag = abs(hx1%tnvec(3)) 
      message = .false.
      ! if ( nzmag .lt. 1.0d0 )  then
         hxArea2D        = triArea2D ( hx1%tvtx%vertexA(1:2), hx1%tvtx%vertexB(1:2), hx1%tvtx%vertexC(1:2) )
         pArea(1)        = triArea2D ( rayvec(1:2),           hx1%tvtx%vertexA(1:2), hx1%tvtx%vertexB(1:2) )
         pArea(2)        = triArea2D ( rayvec(1:2),           hx1%tvtx%vertexB(1:2), hx1%tvtx%vertexC(1:2) )
         pArea(3)        = triArea2D ( rayvec(1:2),           hx1%tvtx%vertexC(1:2), hx1%tvtx%vertexA(1:2) )
         sumArea = pArea(1) + pArea(2) +  pArea(3) 
         if (sumArea .le. hxArea2D *(1.0d0 + tiny))  message = .true. 
         ! 
   end subroutine 
   !*********************************************************************
   !  check_ray_triangle_intersection - to check a point is within a triangle.
   !*********************************************************************
   subroutine check_ray_triangle_intersection  & 
        (rayvec, hx1, nout, verbose, message)
      implicit none
      real(8)           , intent(in)    :: rayvec(3)
      type(stl_tri)     , intent(in)    :: hx1
      integer           , intent(in)    :: nout
      logical           , intent(in)    :: verbose
      logical           , intent(inout) :: message
      real(8)   :: normalVec(3), pArea(3), sumArea, tiny=1.E-2, nzmag
      real(8)   :: pseudoParaAreaSquare, tinyNum=0.0E-3
      nzmag = abs(hx1%tnvec(3)) 
      message = .false.
      if ( nzmag < (1.00d0 - tinyNum) )  then 
         call calc_tri_area (rayvec,hx1%tvtx%vertexA,hx1%tvtx%vertexB,pArea(1))
         call calc_tri_area (rayvec,hx1%tvtx%vertexB,hx1%tvtx%vertexC,pArea(2))
         call calc_tri_area (rayvec,hx1%tvtx%vertexC,hx1%tvtx%vertexA,pArea(3))
         sumArea = pArea(1) + pArea(2) +  pArea(3) 
         sumArea = sumArea * (1.0D0 - tiny*1.0d0)
         if (sumArea .le. hx1%area )  message = .true. 
      end if
   end subroutine check_ray_triangle_intersection
   !*********************************************************************
   ! update_zsurf_top_btm - to update the top and btm position (z-coord)
   !    at a specific grid point (ipx, jpy).
   !*********************************************************************
   subroutine UPDATE_ZSURF_TOP_BTM  & 
        (itri, rayvec, zval_btm, zval_top, id_btm, id_top)
      implicit none
      integer   , intent(in)    :: itri
      integer   , intent(inout) :: id_btm, id_top
      real(8)   , intent(in)    :: rayvec(3)
      real(8)   , intent(inout) :: zval_btm, zval_top
      real(8)                   :: tmpval
      if      (rayvec(3) < zval_btm .and. rayvec(3) > zval_top ) then
         zval_btm  = rayvec(3) 
         zval_top  = rayvec(3) 
         id_top = itri
         id_btm = itri
      else if (rayvec(3) < zval_btm  ) then
         ! zval_top   = zval_btm
         zval_btm  = rayvec(3)
         id_btm = itri
      else if (rayvec(3) > zval_top  ) then
         zval_top   = rayvec(3)
         id_top = itri
      else       ! not important case
      endif
      if (zval_btm > zval_top)  then
         tmpval= zval_top
         zval_top   = zval_btm
         zval_btm   = tmpval
      endif
   end subroutine UPDATE_ZSURF_TOP_BTM
   !*********************************************************************
   ! update_zsurf_top_btm - to update the top and btm position (z-coord)
   !    at a specific grid point (ipx, jpy).
   !*********************************************************************
   subroutine UPDATE_ZSURF_TOP_BTM_v0  & 
        (itri, rayvec, zval_btm, zval_top, id_btm, id_top)
      implicit none
      integer   , intent(in)    :: itri
      integer   , intent(inout) :: id_btm, id_top
      real(8)   , intent(in)    :: rayvec(3)
      real(8)   , intent(inout) :: zval_btm, zval_top
      real(8)                   :: tmpval
      ! write(*,*) "Ha Ha "
      if      (rayvec(3) < zval_btm .and. rayvec(3) > zval_top ) then
         zval_btm  = rayvec(3) 
         zval_top  = rayvec(3) 
         id_top = itri
         id_btm = itri
      else if (rayvec(3) < zval_btm  ) then
         zval_btm  = rayvec(3)
         id_btm = itri
      else if (rayvec(3) > zval_top  ) then
         zval_top   = rayvec(3)
         id_top = itri
      else       ! not important case
      endif
      if (zval_btm > zval_top)  then
         tmpval= zval_top
         zval_top   = zval_btm
         zval_btm   = tmpval
      endif
   end subroutine UPDATE_ZSURF_TOP_BTM_v0
   !*********************************************************************
   ! update_zsurf_top_btm2 - to update the top and btm position (z-coord)
   !    at a specific grid point (ipx, jpy).
   !*********************************************************************
   subroutine UPDATE_ZSURF_TOP_BTM2 & 
        (itri, rayvec, zval_btm, zval_top, id_btm, id_top)
      implicit none
      integer   , intent(in)    :: itri
      integer   , intent(inout) :: id_btm, id_top
      real(8)   , intent(in)    :: rayvec(3)
      real(8)   , intent(inout) :: zval_btm, zval_top
      real(8)                   :: tmpval
      ! write(*,*) "Ha Ha "
      if      (rayvec(3) < zval_btm .and. rayvec(3) > zval_top ) then
         ! first value, go to the bottom
         zval_btm  = rayvec(3) 
         zval_top  = rayvec(3) 
         id_top = itri
         id_btm = itri
      else if (rayvec(3) < zval_btm  ) then
         ! zval_top   = zval_btm
         zval_btm  = rayvec(3)
         id_btm = itri
      else if (rayvec(3) > zval_top  ) then
         zval_top   = rayvec(3)
         id_top = itri
      else       ! not important case
      endif
      ! write(*,*) "Ha Ha 2"
      if (zval_btm > zval_top)  then
         tmpval= zval_top
         zval_top   = zval_btm
         zval_btm   = tmpval
      endif
   end subroutine UPDATE_ZSURF_TOP_BTM2
   !*********************************************************************
   ! calc_a_tri_nvec - to calculate a normal vector of a triangle
   !*********************************************************************
   subroutine CALC_A_TRI_NVEC ( Nparity,hx1 )
      implicit none 
      type(stl_tri), intent(inout) :: hx1
      integer   :: Nparity   
      real (8)  :: nmlvector(3)
      hx1%trvs%vecBwrtA  = hx1%tvtx%vertexB - hx1%tvtx%vertexA 
      hx1%trvs%vecCwrtB  = hx1%tvtx%vertexC - hx1%tvtx%vertexB
      hx1%trvs%vecAwrtC  = hx1%tvtx%vertexA - hx1%tvtx%vertexC
      nmlvector = (hx1%trvs%vecBwrtA .cross. hx1%trvs%vecCwrtB) !  2018-07-21-18-42-20-PM
      hx1%tnvec = uvec(nmlvector) * DBLE(Nparity)
   end subroutine CALC_A_TRI_NVEC
   !*********************************************************************
   ! make_stl_triangles_top
   ! make_stl_triangles_frt_bck
   ! make_stl_triangles_rit_lft
   !*********************************************************************
   subroutine MAKE_STL_TRIANGLES_TOP & 
        (Npx, Npy, Ndx, Ndy, postnX, postnY, z_surf_top, grhx_top, zgrid_shift )
      implicit none
      ! STL triangles on the btm surfaces can be calculated by setting the parity = - parity.
      ! 2018-07-21-11-58-37-AM
      integer, intent(in) :: Npx, Npy, Ndx, Ndy
      real(8), intent(in) :: postnX(Npx), postnY (Npy)
      real(8), intent(in) :: z_surf_top(Npx, Npy), zgrid_shift
      type(stl_tri) :: grhx_top (Ndx*Ndy*2)
      integer :: ipx, jpy, ktri
      integer :: i1, j1, i2,j2, i3,j3
      integer :: ipx_self, jpy_self 
      integer :: ipx_plus, jpy_plus
      integer :: ipx_mnus, jpy_mnus
      ktri = 0
      do ipx = 2, Ndx, 2
         do jpy = 2, Ndy, 2
            ipx_self = ipx ; ipx_plus = ipx + 1 ; ipx_mnus = ipx - 1
            jpy_self = jpy ; jpy_plus = jpy + 1 ; jpy_mnus = jpy - 1
            ! 1. Eeast-East-North
               i1 = ipx_self ; j1 = jpy_self ; 
               i2 = ipx_plus ; j2 = jpy_self ; 
               i3 = ipx_plus ; j3 = jpy_plus
               ktri = ktri + 1
               grhx_top(ktri)%tvtx%vertexA(1) = postnX     (i1)
               grhx_top(ktri)%tvtx%vertexA(2) = postnY     (j1)
               grhx_top(ktri)%tvtx%vertexB(1) = postnX     (i2)
               grhx_top(ktri)%tvtx%vertexB(2) = postnY     (j2)
               grhx_top(ktri)%tvtx%vertexC(1) = postnX     (i3)
               grhx_top(ktri)%tvtx%vertexC(2) = postnY     (j3)
               grhx_top(ktri)%tvtx%vertexA(3) = z_surf_top (i1,j1)
               grhx_top(ktri)%tvtx%vertexB(3) = z_surf_top (i2,j2) 
               grhx_top(ktri)%tvtx%vertexC(3) = z_surf_top (i3,j3)
            ! 2. East-North-North
               i1 = ipx_self ; j1 = jpy_self ; 
               i2 = ipx_plus ; j2 = jpy_plus ; 
               i3 = ipx_self ; j3 = jpy_plus
               ktri = ktri + 1
               grhx_top(ktri)%tvtx%vertexA(1) = postnX     (i1)
               grhx_top(ktri)%tvtx%vertexA(2) = postnY     (j1)
               grhx_top(ktri)%tvtx%vertexB(1) = postnX     (i2)
               grhx_top(ktri)%tvtx%vertexB(2) = postnY     (j2)
               grhx_top(ktri)%tvtx%vertexC(1) = postnX     (i3)
               grhx_top(ktri)%tvtx%vertexC(2) = postnY     (j3)
               grhx_top(ktri)%tvtx%vertexA(3) = z_surf_top (i1,j1)
               grhx_top(ktri)%tvtx%vertexB(3) = z_surf_top (i2,j2) 
               grhx_top(ktri)%tvtx%vertexC(3) = z_surf_top (i3,j3)
            ! 3. West-North-North
               i1 = ipx_self ; j1 = jpy_self ; 
               i2 = ipx_self ; j2 = jpy_plus ; 
               i3 = ipx_mnus ; j3 = jpy_plus
               ktri = ktri + 1
               grhx_top(ktri)%tvtx%vertexA(1) = postnX     (i1)
               grhx_top(ktri)%tvtx%vertexA(2) = postnY     (j1)
               grhx_top(ktri)%tvtx%vertexB(1) = postnX     (i2)
               grhx_top(ktri)%tvtx%vertexB(2) = postnY     (j2)
               grhx_top(ktri)%tvtx%vertexC(1) = postnX     (i3)
               grhx_top(ktri)%tvtx%vertexC(2) = postnY     (j3)
               grhx_top(ktri)%tvtx%vertexA(3) = z_surf_top (i1,j1)
               grhx_top(ktri)%tvtx%vertexB(3) = z_surf_top (i2,j2) 
               grhx_top(ktri)%tvtx%vertexC(3) = z_surf_top (i3,j3)
            ! 4. West-West-North
               i1 = ipx_self ; j1 = jpy_self ; 
               i2 = ipx_mnus ; j2 = jpy_plus ; 
               i3 = ipx_mnus ; j3 = jpy_self
               ktri = ktri + 1
               grhx_top(ktri)%tvtx%vertexA(1) = postnX     (i1)
               grhx_top(ktri)%tvtx%vertexA(2) = postnY     (j1)
               grhx_top(ktri)%tvtx%vertexB(1) = postnX     (i2)
               grhx_top(ktri)%tvtx%vertexB(2) = postnY     (j2)
               grhx_top(ktri)%tvtx%vertexC(1) = postnX     (i3)
               grhx_top(ktri)%tvtx%vertexC(2) = postnY     (j3)
               grhx_top(ktri)%tvtx%vertexA(3) = z_surf_top (i1,j1)
               grhx_top(ktri)%tvtx%vertexB(3) = z_surf_top (i2,j2) 
               grhx_top(ktri)%tvtx%vertexC(3) = z_surf_top (i3,j3)
            ! 5. West-West-South
               i1 = ipx_self ; j1 = jpy_self ; 
               i2 = ipx_mnus ; j2 = jpy_self ; 
               i3 = ipx_mnus ; j3 = jpy_mnus
               ktri = ktri + 1
               grhx_top(ktri)%tvtx%vertexA(1) = postnX     (i1)
               grhx_top(ktri)%tvtx%vertexA(2) = postnY     (j1)
               grhx_top(ktri)%tvtx%vertexB(1) = postnX     (i2)
               grhx_top(ktri)%tvtx%vertexB(2) = postnY     (j2)
               grhx_top(ktri)%tvtx%vertexC(1) = postnX     (i3)
               grhx_top(ktri)%tvtx%vertexC(2) = postnY     (j3)
               grhx_top(ktri)%tvtx%vertexA(3) = z_surf_top (i1,j1)
               grhx_top(ktri)%tvtx%vertexB(3) = z_surf_top (i2,j2) 
               grhx_top(ktri)%tvtx%vertexC(3) = z_surf_top (i3,j3)
            ! 6. West-South-South
               i1 = ipx_self ; j1 = jpy_self ; 
               i2 = ipx_mnus ; j2 = jpy_mnus ; 
               i3 = ipx_self ; j3 = jpy_mnus
               ktri = ktri + 1
               grhx_top(ktri)%tvtx%vertexA(1) = postnX     (i1)
               grhx_top(ktri)%tvtx%vertexA(2) = postnY     (j1)
               grhx_top(ktri)%tvtx%vertexB(1) = postnX     (i2)
               grhx_top(ktri)%tvtx%vertexB(2) = postnY     (j2)
               grhx_top(ktri)%tvtx%vertexC(1) = postnX     (i3)
               grhx_top(ktri)%tvtx%vertexC(2) = postnY     (j3)
               grhx_top(ktri)%tvtx%vertexA(3) = z_surf_top (i1,j1)
               grhx_top(ktri)%tvtx%vertexB(3) = z_surf_top (i2,j2) 
               grhx_top(ktri)%tvtx%vertexC(3) = z_surf_top (i3,j3)
            ! 7. East-South-South
               i1 = ipx_self ; j1 = jpy_self ; 
               i2 = ipx_self ; j2 = jpy_mnus ; 
               i3 = ipx_plus ; j3 = jpy_mnus
               ktri = ktri + 1
               grhx_top(ktri)%tvtx%vertexA(1) = postnX     (i1)
               grhx_top(ktri)%tvtx%vertexA(2) = postnY     (j1)
               grhx_top(ktri)%tvtx%vertexB(1) = postnX     (i2)
               grhx_top(ktri)%tvtx%vertexB(2) = postnY     (j2)
               grhx_top(ktri)%tvtx%vertexC(1) = postnX     (i3)
               grhx_top(ktri)%tvtx%vertexC(2) = postnY     (j3)
               grhx_top(ktri)%tvtx%vertexA(3) = z_surf_top (i1,j1)
               grhx_top(ktri)%tvtx%vertexB(3) = z_surf_top (i2,j2) 
               grhx_top(ktri)%tvtx%vertexC(3) = z_surf_top (i3,j3)
            ! 8. East-East-South
               i1 = ipx_self ; j1 = jpy_self ; 
               i2 = ipx_plus ; j2 = jpy_mnus ; 
               i3 = ipx_plus ; j3 = jpy_self 
               ktri = ktri + 1
               grhx_top(ktri)%tvtx%vertexA(1) = postnX     (i1)
               grhx_top(ktri)%tvtx%vertexA(2) = postnY     (j1)
               grhx_top(ktri)%tvtx%vertexB(1) = postnX     (i2)
               grhx_top(ktri)%tvtx%vertexB(2) = postnY     (j2)
               grhx_top(ktri)%tvtx%vertexC(1) = postnX     (i3)
               grhx_top(ktri)%tvtx%vertexC(2) = postnY     (j3)
               grhx_top(ktri)%tvtx%vertexA(3) = z_surf_top (i1,j1)
               grhx_top(ktri)%tvtx%vertexB(3) = z_surf_top (i2,j2) 
               grhx_top(ktri)%tvtx%vertexC(3) = z_surf_top (i3,j3)
               ! 
         end do
      end do
      if (zgrid_shift .ne. 0.0D0)  then
         grhx_top%tvtx%vertexA(3) = grhx_top%tvtx%vertexA(3) + zgrid_shift
         grhx_top%tvtx%vertexB(3) = grhx_top%tvtx%vertexB(3) + zgrid_shift
         grhx_top%tvtx%vertexC(3) = grhx_top%tvtx%vertexC(3) + zgrid_shift
      end if
   end subroutine MAKE_STL_TRIANGLES_TOP
   !*********************************************************************
   subroutine MAKE_STL_TRIANGLES_FRT_BCK (& 
         jpy_dflt, Npx, Npy, Ndx, Ndy, postnX, postnY, z_surf_top, z_surf_btm, & 
        grhx_top, grhx_btm, grhx_frt)
      implicit none
      integer, intent(in) :: jpy_dflt
      integer, intent(in) :: Npx, Npy, Ndx, Ndy
      real(8), intent(in) :: postnX(Npx), postnY (Npy)
      real(8), intent(in) :: z_surf_top(Npx, Npy),z_surf_btm(Npx, Npy)
      type(stl_tri)     :: grhx_top (Ndx*Ndy*2), grhx_btm(Ndx*Ndy*2), grhx_4sd((Ndx+Ndy)*2*2)
      type(stl_tri)     :: grhx_frt(Ndx*2)
      integer :: i1, j1, i2, j2, i3, j3, i4, j4
      integer :: ipx, jpy, ktri
      integer :: ipx_self, jpy_self 
      integer :: ipx_plus, jpy_plus
      integer :: ipx_mnus, jpy_mnus
      j1 = jpy_dflt
      ktri = 0 
      DO ipx = 2 , Ndx, 2
         ipx_self = ipx ; ipx_plus = ipx + 1 ; ipx_mnus = ipx - 1
         ! triangle 1: b2 b3 t3 =====================
         ktri = ktri +1
         ! b2 point
         i1 = ipx_self ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexA(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexA(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexA(3) = z_surf_btm (i1,j1)
         ! b3 point
         i1 = ipx_plus ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexB(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexB(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexB(3) = z_surf_btm (i1,j1)
         ! t3 point
         i1 = ipx_plus ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexC(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexC(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexC(3) = z_surf_top (i1,j1)
         ! triangle 2: b2 t3 t2 =====================
         ktri = ktri +1
         ! b2 point
         i1 = ipx_self ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexA(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexA(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexA(3) = z_surf_btm (i1,j1)
         ! t3 point
         i1 = ipx_plus ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexB(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexB(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexB(3) = z_surf_top (i1,j1)
         ! t2 point
         i1 = ipx_self ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexC(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexC(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexC(3) = z_surf_top (i1,j1)
         ! triangle 3: b2 t2 t1 =====================
         ktri = ktri +1
         ! b2 point
         i1 = ipx_self ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexA(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexA(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexA(3) = z_surf_btm (i1,j1)
         ! t2 point
         i1 = ipx_self ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexB(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexB(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexB(3) = z_surf_top (i1,j1)
         ! t1 point
         i1 = ipx_mnus ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexC(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexC(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexC(3) = z_surf_top (i1,j1)
         ! triangle 4: b2 t1 b1 =====================
         ktri = ktri +1
         ! b2 point
         i1 = ipx_self ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexA(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexA(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexA(3) = z_surf_btm (i1,j1)
         ! t1 point
         i1 = ipx_mnus ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexB(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexB(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexB(3) = z_surf_top (i1,j1)
         ! b1 point
         i1 = ipx_mnus ; j1 = jpy_dflt
            grhx_frt(ktri)%tvtx%vertexC(1) = postnX     (i1)
            grhx_frt(ktri)%tvtx%vertexC(2) = postnY     (j1)
            grhx_frt(ktri)%tvtx%vertexC(3) = z_surf_btm (i1,j1)
         ! 
      END DO
      return
   end subroutine MAKE_STL_TRIANGLES_FRT_BCK
   !*********************************************************************
   subroutine MAKE_STL_TRIANGLES_RIT_LFT ( & 
        ipx_dflt,Npx, Npy, Ndx, Ndy, postnX, postnY, z_surf_top, z_surf_btm, & 
        grhx_top, grhx_btm, grhx_rit)
      implicit none
      integer, intent(in) :: ipx_dflt
      integer, intent(in) :: Npx, Npy, Ndx, Ndy
      real(8), intent(in) :: postnX(Npx), postnY (Npy)
      real(8), intent(in) :: z_surf_top(Npx, Npy),z_surf_btm(Npx, Npy)
      type(stl_tri)     :: grhx_top (Ndx*Ndy*2), grhx_btm(Ndx*Ndy*2)
      type(stl_tri)     :: grhx_rit(Ndy*2)
      integer :: i1, j1, i2, j2, i3, j3, i4, j4
      integer :: ipx, jpy, ktri
      integer :: ipx_self, jpy_self 
      integer :: ipx_plus, jpy_plus
      integer :: ipx_mnus, jpy_mnus
      i1 = ipx_dflt
      ktri = 0 
      DO jpy = 2 , Ndy, 2
         jpy_self = jpy ; jpy_plus = jpy + 1 ; jpy_mnus = jpy - 1
         ! triangle 1: b2 b3 t3 =====================
         ktri = ktri +1
         ! b2 point
         i1 = ipx_dflt ; j1 = jpy_self
            grhx_rit(ktri)%tvtx%vertexA(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexA(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexA(3) = z_surf_btm (i1,j1)
         ! b3 point
         i1 = ipx_dflt ; j1 = jpy_plus
            grhx_rit(ktri)%tvtx%vertexB(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexB(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexB(3) = z_surf_btm (i1,j1)
         ! t3 point
         i1 = ipx_dflt ; j1 = jpy_plus
            grhx_rit(ktri)%tvtx%vertexC(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexC(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexC(3) = z_surf_top (i1,j1)
         ! triangle 2: b2 t3 t2 =====================
         ktri = ktri +1
         ! b2 point
         i1 = ipx_dflt ; j1 = jpy_self
            grhx_rit(ktri)%tvtx%vertexA(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexA(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexA(3) = z_surf_btm (i1,j1)
         ! t3 point
         i1 = ipx_dflt ; j1 = jpy_plus
            grhx_rit(ktri)%tvtx%vertexB(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexB(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexB(3) = z_surf_top (i1,j1)
         ! t2 point
         i1 = ipx_dflt ; j1 = jpy_self
            grhx_rit(ktri)%tvtx%vertexC(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexC(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexC(3) = z_surf_top (i1,j1)
         ! triangle 3: b2 t2 t1 =====================
         ktri = ktri +1
         ! b2 point
         i1 = ipx_dflt ; j1 = jpy_self
            grhx_rit(ktri)%tvtx%vertexA(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexA(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexA(3) = z_surf_btm (i1,j1)
         ! t2 point
         i1 = ipx_dflt ; j1 = jpy_self
            grhx_rit(ktri)%tvtx%vertexB(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexB(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexB(3) = z_surf_top (i1,j1)
         ! t1 point
         i1 = ipx_dflt ; j1 = jpy_mnus
            grhx_rit(ktri)%tvtx%vertexC(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexC(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexC(3) = z_surf_top (i1,j1)
         ! triangle 4: b2 t1 b1 =====================
         ktri = ktri +1
         ! b2 point
         i1 = ipx_dflt ; j1 = jpy_self
            grhx_rit(ktri)%tvtx%vertexA(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexA(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexA(3) = z_surf_btm (i1,j1)
         ! t1 point
         i1 = ipx_dflt ; j1 = jpy_mnus
            grhx_rit(ktri)%tvtx%vertexB(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexB(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexB(3) = z_surf_top (i1,j1)
         ! b1 point
         i1 = ipx_dflt ; j1 = jpy_mnus
            grhx_rit(ktri)%tvtx%vertexC(1) = postnX     (i1)
            grhx_rit(ktri)%tvtx%vertexC(2) = postnY     (j1)
            grhx_rit(ktri)%tvtx%vertexC(3) = z_surf_btm (i1,j1)
      END DO
   end subroutine MAKE_STL_TRIANGLES_RIT_LFT
   !*********************************************************************
   !*********************************************************************
   ! smooth_blank_nodes
   !*********************************************************************
   subroutine smooth_blank_nodes ( &
         Npx, Npy, z_surf_top, ncount_zero_top, list_zero_top, &
         z_surf_btm, ncount_zero_btm, list_zero_btm  )
      implicit none
      integer, intent(in) :: Npx, Npy
      integer, intent(in) :: ncount_zero_top
      integer, intent(in) :: ncount_zero_btm
      type(stl_list),   intent(in)      :: list_zero_top(ncount_zero_top)
      type(stl_list),   intent(in)      :: list_zero_btm(ncount_zero_btm)
      real(8),          intent(inout)   :: z_surf_btm(Npx, Npy)
      real(8),          intent(inout)   :: z_surf_top(Npx, Npy)
      ! 
      integer   :: ipx, ipx_plus, ipx_mnus
      integer   :: jpy, jpy_plus, jpy_mnus
      integer   :: itri, nNbor_top, nNbor_btm 
      integer   :: iMC
      real(8)   :: val_surf_top, val_surf_btm
      ! print *, ncount_zero_top,  ncount_zero_btm
      do itri = 1, ncount_zero_top
         ipx = list_zero_top(itri)%ix
         jpy = list_zero_top(itri)%jy
         z_surf_top(ipx, jpy) = 0.0
      end do
      do itri = 1, ncount_zero_btm
         ipx = list_zero_btm(itri)%ix
         jpy = list_zero_btm(itri)%jy
         z_surf_btm(ipx, jpy) = 0.0
      end do
      iMC_loop: do iMC = 1, ncount_zero_top * 1000
         itir_top_loop: do itri = 1, ncount_zero_top
            ipx            = list_zero_top(itri)%ix
            jpy            = list_zero_top(itri)%jy
            ipx_plus       = ipx + 1
            ipx_mnus       = ipx - 1
            jpy_plus       = jpy + 1
            jpy_mnus       = jpy - 1
            nNbor_top      = 0
            ! update Top
            val_surf_top = 0
            if (ipx_plus .le. Npx) then
               val_surf_top = val_surf_top + z_surf_top(ipx_plus, jpy) 
               nNbor_top = nNbor_top + 1
            end if
            if (ipx_mnus .ge. 1) then
               val_surf_top = val_surf_top + z_surf_top(ipx_mnus, jpy) 
               nNbor_top = nNbor_top + 1
            end if
            if (jpy_plus .le. Npy) then
               val_surf_top = val_surf_top + z_surf_top(ipx, jpy_plus) 
               nNbor_top = nNbor_top + 1
            end if
            if (jpy_mnus .ge. 1) then
               val_surf_top = val_surf_top + z_surf_top(ipx, jpy_mnus) 
               nNbor_top = nNbor_top + 1
            end if
            z_surf_top(ipx, jpy)  =  val_surf_top / dble(nNbor_top)
         end do itir_top_loop
         itir_btm_loop: do itri = 1, ncount_zero_btm
            ipx            = list_zero_btm(itri)%ix
            jpy            = list_zero_btm(itri)%jy
            ipx_plus       = ipx + 1
            ipx_mnus       = ipx - 1
            jpy_plus       = jpy + 1
            jpy_mnus       = jpy - 1
            nNbor_btm      = 0
            ! update Btm
            val_surf_btm = 0
            if (ipx_plus .le. Npx) then
               val_surf_btm = val_surf_btm + z_surf_btm(ipx_plus, jpy) 
               nNbor_btm = nNbor_btm + 1
            end if
            if (ipx_mnus .ge. 1) then
               val_surf_btm = val_surf_btm + z_surf_btm(ipx_mnus, jpy) 
               nNbor_btm = nNbor_btm + 1
            end if
            if (jpy_plus .le. Npy) then
               val_surf_btm = val_surf_btm + z_surf_btm(ipx, jpy_plus) 
               nNbor_btm = nNbor_btm + 1
            end if
            if (jpy_mnus .ge. 1) then
               val_surf_btm = val_surf_btm + z_surf_btm(ipx, jpy_mnus) 
               nNbor_btm = nNbor_btm + 1
            end if
            z_surf_btm(ipx, jpy)  =  val_surf_btm / dble(nNbor_btm)
         end do itir_btm_loop
      end do iMC_loop
   end subroutine smooth_blank_nodes
!*******************************************************************
!
!                            PRINT/WRITE ROUTINES
! 
!*******************************************************************
   !************************************************************************
   !  print_info of a triangle with its index, cid is calculated by default.
   !************************************************************************
   subroutine print_a_triangle (itri,hx1,nout,verbose)     
      implicit none
      type(stl_tri), intent(in) :: hx1
      integer                   :: itri,cid,nout
      logical                   :: verbose
      ! cid            = check_uvec_dvec (hx1%tnvec)
      if (verbose .eqv. .true.)  then 
         write(nout,FMT_1C1I) "nid ", itri
         write(nout,FMT_1C1I) "cid ", hx1%cid
         write(nout,FMT_1C4RN)"nvec", hx1%tnvec(:)
         write(nout,FMT_1C4RN)"vrtA", hx1%tvtx%vertexA
         write(nout,FMT_1C4RN)"vrtB", hx1%tvtx%vertexB
         write(nout,FMT_1C4RN)"vrtC", hx1%tvtx%vertexC
         write(nout,FMT_1C4RN)"tmid", hx1%tmid
         write(nout,FMT_1C4RN)"vB-A", hx1%trvs%vecBwrtA
         write(nout,FMT_1C4RN)"vC-B", hx1%trvs%vecCwrtB
         write(nout,FMT_1C4RN)"vA-C", hx1%trvs%vecAwrtC
         write(nout,FMT_1C4RN)"area", hx1%area
      end if
      ! 
   end subroutine print_a_triangle
   !****************************************************************
   !  print_info - print out detailed interpolation info 
   !****************************************************************
   subroutine print_info (itri,hx1,mat44a,mat44b,det44a,det44b,zpostn,nout,verbose)
      ! This is to print out information of 4x4 matrix and interpolated value of zmid.
      ! 2018-07-19-08-13-27-AM
      implicit none
      integer, intent(in)       :: itri 
      integer                   :: cid,nout
      type(stl_tri), intent(in) :: hx1
      real(8) :: mat44a(4,4),mat44b(4,4),det44a,det44b,zpostn
      logical                   :: verbose
      cid            = check_uvec_dvec (hx1%tnvec)
      ! write(*,*) "===========print_info"
      call  print_a_triangle (itri,hx1,nout,verbose)     
      return
   end subroutine print_info
   !****************************************************************
   !  write_ascii_stl - to write an ascii stl file 
   !****************************************************************
   subroutine write_ascii_stl (outputfile,Nx,Ny,hx,title)
      implicit none 
      integer,                  intent(in)      :: Nx,Ny
      type(stl_tri),            intent(inout)   :: hx(Nx,Ny,2)
      character(len=128),       intent(inout)   :: outputfile, title
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat, lin=10, i,j ,k
      OPEN(lin, file=outputfile, status='replace')
         WRITE(lin,"('solid',2X,A128)") title 
         write(*,*) title
         ix_loop: do i = 1, Nx-1
            jy_loop: do j = 1, Nx-1
               k_loop: do k = 1,2
               enddo k_loop
            enddo jy_loop
         enddo ix_loop
         WRITE(lin,"('endsolid',2X,A128)") title 
      CLOSE(lin)
   end subroutine write_ascii_stl
   !****************************************************************
   !  write_ascii_stl - to write an ascii stl file 
   !****************************************************************
   subroutine write_ascii_stl_each (outputfile,lin,hx1)
      implicit none 
      type(stl_tri),            intent(inout)   :: hx1 
      character(len=128),       intent(inout)   :: outputfile
      integer   , intent(in) :: lin
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat,  i,j 
                  WRITE(lin,"(2x,'facet normal',1X,3(1x,ES16.8))")  hx1%tnvec
                  WRITE(lin,"(4x,'outer loop')")
                  WRITE(lin,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexA
                  WRITE(lin,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexB
                  WRITE(lin,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexC
                  WRITE(lin,"(4x,'endloop')")
                  WRITE(lin,"(4x,'endfacet')")
   end subroutine 
   !****************************************************************
   ! write_blank_nodes
   !****************************************************************
   subroutine write_blank_nodes ( Npx, Npy, Max_Z, Min_Z, & 
        nchkTop, z_surf_top, chkfile_top00, & 
        nchkBtm, z_surf_btm, chkfile_btm00, &
        ncount_zero_top, ncount_zero_btm)
      implicit none
      character(len=128), intent(in) :: chkfile_top00, chkfile_btm00
      integer, intent(in) :: Npx, Npy, nchkTop, nchkBtm
      real(8), intent(in) :: z_surf_btm(Npx, Npy), z_surf_top(Npx, Npy)
      real(8), intent(in) :: Min_Z, Max_Z 
      integer, intent(out)::  ncount_zero_top, ncount_zero_btm
      integer :: ipx, jpy
      open(nchkTop,  file=chkfile_top00,  status="replace")
      open(nchkBtm,  file=chkfile_btm00,  status="replace")
      ncount_zero_top = 0 ; ncount_zero_btm = 0
      ipx_loop: do ipx = 1, Npx
         jpy_loop: do jpy = 1, Npy
            if ( z_surf_top(ipx,jpy) == 0) then
               ncount_zero_top = ncount_zero_top + 1
               write(nchkTop,*) ipx, jpy!, id_tri_btm(ipx,jpy)
            end if
            if ( z_surf_btm(ipx,jpy) > Max_Z ) then
               ncount_zero_btm = ncount_zero_btm + 1
               write(nchkBtm,*) ipx, jpy! , id_tri_top(ipx,jpy)
            end if
         end do jpy_loop
      end do ipx_loop
      close(nchkTop)
      close(nchkBtm)
   end subroutine write_blank_nodes
   !****************************************************************
   !****************************************************************
   subroutine write_a_ascii_stl_tri (nfile,hx1)
      implicit none 
      type(stl_tri),            intent(in)   :: hx1 
      integer   , intent(in) :: nfile
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat,  i,j 
      WRITE(nfile,"(2x,'facet normal',1X,3(1x,ES16.8))")  hx1%tnvec
      WRITE(nfile,"(4x,'outer loop')")
      WRITE(nfile,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexA
      WRITE(nfile,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexB
      WRITE(nfile,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexC
      WRITE(nfile,"(4x,'endloop')")
      WRITE(nfile,"(4x,'endfacet')")
   end subroutine write_a_ascii_stl_tri
! 
   subroutine write_stl_header_footer (stlfid, stlTitle, HoF)
      implicit none
      integer,             intent(in)      :: stlfid
      character*128,       intent(in)      :: stlTitle
      character*1,         intent(in)      :: HoF  ! head or foot
      if (HoF == "H" .or. HoF == "h") then
         WRITE(stlfid,"(  'solid',5X,A12)") trim(stlTitle)
      elseif (HoF == "F" .or. HoF == "f") then
         WRITE(stlfid,"('endsolid',2X,A12)") trim(stlTitle)
      end if
   end subroutine write_stl_header_footer
!
!*******************************************************************
! 
!                            READING ROUTINS
! 
!*******************************************************************
   subroutine read_stl_header_footer (stlfid, stlTitle, HoF ) 
      implicit none
      integer,          intent(in)      :: stlfid
      character*128,    intent(out)     :: stlTitle
      character*1,      intent(in)      :: HoF  ! head or foot
      character(len=128)                :: Csolid, Cendsolid
      if (HoF == "H" .or. HoF == "h") then
         READ(stlfid,*) Csolid,  stlTitle
      elseif (HoF == "F" .or. HoF == "f") then
         READ(stlfid,*) Cendsolid,  stlTitle
      end if
   end subroutine read_stl_header_footer
   !****************************************************************
   subroutine read_a_ascii_stl_tri (nfile,hx1)
      implicit none 
      type(stl_tri), intent(inout)   :: hx1 
      integer , intent(in) :: nfile
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat,  i,j 
      read(nfile,*) Cfacet, Cnormal ,hx1%tnvec
      read(nfile,*) Couter, Cloop   
      read(nfile,*) Cvertex, hx1%tvtx%vertexA
      read(nfile,*) Cvertex, hx1%tvtx%vertexB
      read(nfile,*) Cvertex, hx1%tvtx%vertexC
      read(nfile,*) Cendloop
      read(nfile,*) Cendfacet
   end subroutine read_a_ascii_stl_tri
   !****************************************************************
   !  read_ascii_stlfile - to read an ascii stl file 
   !****************************************************************
   subroutine read_ascii_stl_file (inputfile,numrec,hx,title,verbose)
      ! This routine is idential to read_ascii_stl except verbose option. 
      integer,                  intent(in)      :: numrec
      type(stl_tri),            intent(inout)   :: hx(numrec)
      logical,                  intent(in)      :: verbose
      character(len=128),       intent(inout)   :: inputfile, title
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat, lin=10
      character(len=128)   :: screenMesssage = & 
        "=== The title of the stl file is (from subroutine read_ascii_stl): "
      OPEN(lin, file=inputfile, status='old')
         READ(lin,"(A128)",iostat=eastat) title 
         if (verbose) write(*,*)  trim(screenMesssage)
         if (verbose) write(*,"(5X,A64)")  adjustl(title)
         do itri = 1,numrec
            read(lin,*) Cfacet, Cnormal ,hx(itri)%tnvec
            read(lin,*) Couter, Cloop   
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexA
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexB
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexC
            read(lin,*) Cendloop
            read(lin,*) Cendfacet
            ! write(*,*)  itri
         enddo 
      CLOSE(lin)
   end subroutine read_ascii_stl_file
   !****************************************************************
   !  read_ascii_stlfile_v2 - to read an ascii stl file 
   !****************************************************************
   subroutine read_ascii_stl_file_v2 (inputfile,numrec,hx,title,verbose)
      ! This routine is idential to read_ascii_stl except verbose option. 
      character(*),             intent(inout)   :: inputfile, title
      integer,                  intent(in)      :: numrec
      type(stl_tri),            intent(inout)   :: hx(numrec)
      logical,                  intent(in)      :: verbose
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat, lin=10
      character(len=128)   :: screenMesssage = & 
        "=== The title of the stl file is (from subroutine read_ascii_stl): "
      OPEN(lin, file=inputfile, status='old')
         READ(lin,"(A128)",iostat=eastat) title 
         if (verbose) write(*,*)  trim(screenMesssage)
         if (verbose) write(*,"(5X,A64)")  adjustl(title)
         do itri = 1,numrec
            read(lin,*) Cfacet, Cnormal ,hx(itri)%tnvec
            read(lin,*) Couter, Cloop   
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexA
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexB
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexC
            read(lin,*) Cendloop
            read(lin,*) Cendfacet
            ! write(*,*)  itri
         enddo 
      CLOSE(lin)
   end subroutine read_ascii_stl_file_v2
   
   !
   
   ! =========== devel starts   2018-08-07-15-35-10-PM
   ! 
   subroutine read_ascii_stl_file_ptr (inputfile,numrec,hx_ptr,title,verbose)
      ! This routine is idential to read_ascii_stl except verbose option. 
      integer,                  intent(in)      :: numrec
      type(stl_tri), pointer,   intent(inout)   :: hx_ptr(:)
      logical,                  intent(in)      :: verbose
      character(len=128),       intent(inout)   :: inputfile, title
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat, lin=10
      character(len=128)   :: screenMesssage = & 
        "=== The title of the stl file is (from subroutine read_ascii_stl): "
      OPEN(lin, file=inputfile, status='old') 
         READ(lin,"(A128)",iostat=eastat) title 
         if (verbose) write(*,*)  trim(screenMesssage)
         if (verbose) write(*,"(5X,A64)")  adjustl(title)
         do itri = 1,numrec
            read(lin,*) Cfacet, Cnormal ,hx_ptr(itri)%tnvec
            read(lin,*) Couter, Cloop   
            read(lin,*) Cvertex, hx_ptr(itri)%tvtx%vertexA
            read(lin,*) Cvertex, hx_ptr(itri)%tvtx%vertexB
            read(lin,*) Cvertex, hx_ptr(itri)%tvtx%vertexC
            read(lin,*) Cendloop
            read(lin,*) Cendfacet
            ! write(*,*)  itri
         enddo 
      CLOSE(lin)
   end subroutine read_ascii_stl_file_ptr
   ! 
   ! =========== devel ends     2018-08-07-15-35-10-PM
   ! 
   !****************************************************************
   !  read_ascii_stl - to read an ascii stl file 
   !****************************************************************
   subroutine read_ascii_stl (inputfile,numrec,hx,title)
      integer,                  intent(in)      :: numrec
      type(stl_tri),            intent(inout)   :: hx(numrec)
      character(len=128),       intent(inout)   :: inputfile, title
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat, lin=10
      character(len=128)   :: screenMesssage = & 
        "=== The title of the stl file is (from subroutine read_ascii_stl): "
      OPEN(lin, file=inputfile, status='old')
         READ(lin,"(A128)",iostat=eastat) title 
         write(*,*)  trim(screenMesssage)
         write(*,"(5X,A64)")  trim(adjustl(title))
         do itri = 1,numrec
            hx(itri)%tid = itri
            read(lin,*) Cfacet, Cnormal ,hx(itri)%tnvec
            read(lin,*) Couter, Cloop   
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexA
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexB
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexC
            read(lin,*) Cendloop
            read(lin,*) Cendfacet
         enddo 
      CLOSE(lin)
   end subroutine read_ascii_stl
   !****************************************************************
   !  read_ascii_stl - to read an ascii stl file 
   !****************************************************************
   subroutine read_ascii_stl_w_nborlist &
        (inputfile,numrec,hx,Npx,Npy,dx,dy,Min_X,Min_Y,gx,gy,nnborMax,title,nnlist)
                                        ! (input_fileSTL, numtri,hx,Npx,Npy,dx,dy,title,nnlist) ;!
      integer,                  intent(in)      :: numrec, Npx,Npy,nnborMax
      real(8),                  intent(in)      :: dx, dy, Min_X, Min_Y,gx, gy
      integer,                  intent(inout)   :: nnlist(Npx,Npy,0:nnborMax)
      type(stl_tri),            intent(inout)   :: hx(numrec)
      character(len=128),       intent(inout)   :: inputfile, title
      real (8)          :: xref, yref , dpx, dpy
      character(len=4)  :: Cloop   
      character(len=4)  :: Cfacet  
      character(len=5)  :: Couter  
      character(len=6)  :: Cnormal , Cvertex
      character(len=7)  :: Cendloop
      character(len=8)  :: Cendfacet
      integer           :: ip, jp , iq, jq, ir, jr
      integer           :: numline, cal_num_of_file_lines ,eastat, lin=10
      character(len=128)   :: & 
        screenMesssage="=== The title of the stl file is (from subroutine read_ascii_stl): "
      xref = Min_X + gx - 0.5d0* dx
      yref = Min_Y + gy - 0.5d0* dy
      OPEN(lin, file=inputfile, status='old')
         READ(lin,"(A128)",iostat=eastat) title 
         write(*,*)  trim(screenMesssage)
         write(*,"(5X,A64)")  adjustl(title)
         do itri = 1,numrec
            read(lin,*) Cfacet, Cnormal ,hx(itri)%tnvec
            read(lin,*) Couter, Cloop   
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexA
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexB
            read(lin,*) Cvertex, hx(itri)%tvtx%vertexC
            read(lin,*) Cendloop
            read(lin,*) Cendfacet
            call analyze_a_triangle (hx(itri))  ! to calculate the midpoint and the area
            ! ip = INT ((hx(itri)%tmid(1)-Min_X-marginX) / dx ) + 1 
            ! jp = INT ((hx(itri)%tmid(2)-Min_Y-marginY) / dy ) + 1
            dpx = (hx(itri)%tmid(1)-xref)
            dpy = (hx(itri)%tmid(2)-yref)
            if (dpx < 0.0d0 .or. dpy < 0.0d0) write(*,*) "dpx, dpy", dpx, dpy
            ip = INT ( dpx / dx ) + 1 
            jp = INT ( dpy / dy ) + 1            
            nnlist(ip,jp,0) = nnlist(ip,jp,0) + 1
            nnlist(ip,jp,nnlist(ip,jp,0)) = itri
         enddo 
      CLOSE(lin)
   end subroutine read_ascii_stl_w_nborlist
   !****************************************************************
   ! read_blank_nodes
   !****************************************************************
   subroutine read_blank_nodes ( &
        Max_Z, chkfile_top00, nchkTop, ncount_zero_top, list_zero_top,    &
        Min_Z, chkfile_btm00, nchkBtm, ncount_zero_btm, list_zero_btm )
      implicit none
      integer, intent(in) :: nchkTop, nchkBtm
      integer, intent(in) :: ncount_zero_top
      integer, intent(in) :: ncount_zero_btm
      type(stl_list), intent(inout) :: list_zero_top(ncount_zero_top)
      type(stl_list), intent(inout) :: list_zero_btm(ncount_zero_btm)
      ! real(8), intent(in) :: z_surf_btm(Npx, Npy), z_surf_top(Npx, Npy)
      real(8), intent(in) :: Min_Z, Max_Z 
      character(len=128), intent(in) :: chkfile_top00, chkfile_btm00
      integer :: ipx, jpy, itri
      open(nchkTop,  file=chkfile_top00,  status="old")
      open(nchkBtm,  file=chkfile_btm00,  status="old")
         do itri = 1, ncount_zero_top
            read(nchkTop,*) list_zero_top(itri)%ix, list_zero_top(itri)%jy
         end do
         do itri = 1, ncount_zero_btm
            read(nchkBtm,*) list_zero_btm(itri)%ix, list_zero_btm(itri)%jy
         end do
      close(nchkTop)
      close(nchkBtm)
      return
   end subroutine read_blank_nodes
   !****************************************************************
   ! 
   !****************************************************************
   ! subroutine make_all_facet_zone_ids (nout, cfile, nfacets, hx , ivar, rvar, msg)
   !    implicit none
   !    character(len=128)        :: cfile
   !    integer, intent(in)       :: nout, nfacets
   !    type(stl_tri), intent(inout)   :: hx(nfacets)
   !    integer :: ivar
   !    real*8  :: rvar
   !    logical :: msg
   !    write(nout,*)    
   !    return
   ! end subroutine make_all_facet_zone_ids
   ! 
   !**********************************************************************
   !  make_a_facet_zone_id - to make zoenID of three vertices of a facet.
   !**********************************************************************
   subroutine make_a_facet_zone_id ( hx1, NRx, NRy, Rx, Ry, nout, cfile, ivar, rvar, lmsg) 
      implicit none
      character(len=128), optional        :: cfile
      integer, optional :: nout
      integer, optional :: ivar
      real*8 , optional :: rvar
      logical, optional :: lmsg
      type(stl_tri), intent(inout)   :: hx1 
      integer :: NRx, NRy
      real*8  :: Rx,  Ry  
      hx1%tvtx%zoneidA(1) = INT ( hx1%tvtx%vertexA(1) / Rx ) + 1
      hx1%tvtx%zoneidA(2) = INT ( hx1%tvtx%vertexA(2) / Ry ) + 1 
      hx1%tvtx%zoneidA(0) = Nry*( hx1%tvtx%zoneidA(2) -  1 ) + &
                                  hx1%tvtx%zoneidA(1)
      hx1%tvtx%zoneidB(1) = INT ( hx1%tvtx%vertexB(1) / Rx ) + 1
      hx1%tvtx%zoneidB(2) = INT ( hx1%tvtx%vertexB(2) / Ry ) + 1 
      hx1%tvtx%zoneidB(0) = Nry*( hx1%tvtx%zoneidB(2) -  1 ) + &
                                  hx1%tvtx%zoneidB(1)
      hx1%tvtx%zoneidC(1) = INT ( hx1%tvtx%vertexC(1) / Rx ) + 1
      hx1%tvtx%zoneidC(2) = INT ( hx1%tvtx%vertexC(2) / Ry ) + 1 
      hx1%tvtx%zoneidC(0) = Nry*( hx1%tvtx%zoneidC(2) -  1 ) + &
                                  hx1%tvtx%zoneidC(1)
      hx1%tvtx%zoneidM(1) = INT ( hx1%tvtx%vertexM(1) / Rx ) + 1
      hx1%tvtx%zoneidM(2) = INT ( hx1%tvtx%vertexM(2) / Ry ) + 1 
      hx1%tvtx%zoneidM(0) = Nry*( hx1%tvtx%zoneidM(2) -  1 ) + &
                                  hx1%tvtx%zoneidM(1)
      return
   end subroutine make_a_facet_zone_id
   !************************************************************************
   !  check_inr1_otr2 - spliting when ONE inner and two outer vertices exit.
   !************************************************************************
   logical function check_inr1_otr2 (hx1,postnCutSurf,nmvecCutSurf,signVec) result (lmsg)
      implicit none
      type(stl_tri), intent(inout)   :: hx1 
      real*8, intent(in)  :: postnCutSurf(3), nmvecCutSurf(3)
      real*8, intent(out) :: signVec(3)
      real*8 :: nmlvect (3), crossPt (3), posvect (3), sign_LR
      lmsg = .false.
      nmlvect =  nmvecCutSurf
      crossPt = postnCutSurf
      posvect = hx1%tvtx%vertexA 
      lmsg = check_vertex_side_sign ( nmlvect, crossPt, posvect, sign_LR) 
      signVec(1) =  sign_LR
      posvect = hx1%tvtx%vertexB
      lmsg = check_vertex_side_sign ( nmlvect, crossPt, posvect, sign_LR) 
      signVec(2) =  sign_LR
      posvect = hx1%tvtx%vertexC
      lmsg = check_vertex_side_sign ( nmlvect, crossPt, posvect, sign_LR) 
      signVec(3) =  sign_LR
   end function check_inr1_otr2
   ! 
   !**********************************************************************
   !  check_lft1_rit2 - spliting when ONE vertices on the left
   !**********************************************************************
   logical function check_lft1_rit2 (hx1,postnCutSurf) result (lmsg)
      implicit none
      type(stl_tri), intent(inout)   :: hx1 
      real*8, intent(in)  :: postnCutSurf(3)
      real*8 :: tmpMinimum 
      lmsg = .false.
      ! A|P|BC
      tmpMinimum = min ( hx1%tvtx%vertexB(1), hx1%tvtx%vertexC(1) )
      if      ( hx1%tvtx%vertexA(1) < tmpMinimum        ) then
         if   ( hx1%tvtx%vertexA(1) < postnCutSurf(1) .and. &
                postnCutSurf(1)     < tmpMinimum        ) lmsg = .true. 
      end if
      ! B|P|CA
      tmpMinimum = min ( hx1%tvtx%vertexC(1),hx1%tvtx%vertexA(1) )
      if      ( hx1%tvtx%vertexB(1) < tmpMinimum        ) then
         if   ( hx1%tvtx%vertexB(1) < postnCutSurf(1) .and. &
                postnCutSurf(1)     < tmpMinimum        ) lmsg = .true. 
      end if
      ! C|P|AB
      tmpMinimum = min ( hx1%tvtx%vertexA(1),hx1%tvtx%vertexB(1) )
      if      ( hx1%tvtx%vertexC(1) < tmpMinimum        ) then
         if   ( hx1%tvtx%vertexC(1) < postnCutSurf(1) .and. &
                postnCutSurf(1)     < tmpMinimum        ) lmsg = .true. 
      end if
      return
   end function check_lft1_rit2
   ! 
   !**********************************************************************
   !  check_btm1_top2 - spliting when ONE vertices on the btm 
   !**********************************************************************
   logical function check_btm1_top2 (hx1,postnCutSurf) result (lmsg)
      implicit none
      type(stl_tri), intent(inout)   :: hx1 
      real*8, intent(in)  :: postnCutSurf(3)
      real*8 :: tmpMinimum 
      lmsg = .false.
      ! A|P|BC
      tmpMinimum = min ( hx1%tvtx%vertexB(2), hx1%tvtx%vertexC(2) )
      if      ( hx1%tvtx%vertexA(2) < tmpMinimum        ) then
         if   ( hx1%tvtx%vertexA(2) < postnCutSurf(2) .and. &
                postnCutSurf(2)     < tmpMinimum        ) lmsg = .true. 
      end if
      ! B|P|CA
      tmpMinimum = min ( hx1%tvtx%vertexC(2),hx1%tvtx%vertexA(2) )
      if      ( hx1%tvtx%vertexB(2) < tmpMinimum        ) then
         if   ( hx1%tvtx%vertexB(2) < postnCutSurf(2) .and. &
                postnCutSurf(2)     < tmpMinimum        ) lmsg = .true. 
      end if
      ! C|P|AB
      tmpMinimum = min ( hx1%tvtx%vertexA(2),hx1%tvtx%vertexB(2) )
      if      ( hx1%tvtx%vertexC(2) < tmpMinimum        ) then
         if   ( hx1%tvtx%vertexC(2) < postnCutSurf(2) .and. &
                postnCutSurf(2)     < tmpMinimum        ) lmsg = .true. 
      end if
      return
   end function check_btm1_top2
   !****************************************************************
   ! check_lft2_rit1 - spliting when TWO vertices on the 
   !****************************************************************
   logical function check_lft2_rit1 (hx1,postnCutSurf) result (lmsg)
      implicit none
      type(stl_tri), intent(inout)   :: hx1 
      ! character(len=128) :: cfile
      ! integer, intent(in) :: nout              
      ! integer :: ivar
      real*8, intent(in)  :: postnCutSurf(3)
      real*8 :: tmpMaximum
      ! logical :: lmsg
      lmsg = .false.
      !  C  A | P | B
      tmpMaximum = max ( hx1%tvtx%vertexC(1),hx1%tvtx%vertexA(1) )
      if      ( tmpMaximum      < hx1%tvtx%vertexB(1)   ) then
         if   ( postnCutSurf(1) < hx1%tvtx%vertexB(1)   .and. &
                tmpMaximum      < postnCutSurf(1)       ) lmsg = .true. 
      end if
      !  B  C  | P | A
      tmpMaximum = max ( hx1%tvtx%vertexB(1), hx1%tvtx%vertexC(1) )
      if      ( tmpMaximum      < hx1%tvtx%vertexA(1)   ) then
         if   ( postnCutSurf(1) < hx1%tvtx%vertexA(1)   .and. &
                tmpMaximum      < postnCutSurf(1)       ) lmsg = .true. 
      end if
      !  A  B  | P | C
      tmpMaximum = max ( hx1%tvtx%vertexA(1),hx1%tvtx%vertexB(1) )
      if      ( tmpMaximum      < hx1%tvtx%vertexC(1)   ) then
         if   ( postnCutSurf(1) < hx1%tvtx%vertexC(1)   .and. &
                tmpMaximum      < postnCutSurf(1)       ) lmsg = .true. 
      end if
      return
   end function check_lft2_rit1
   ! 
   !****************************************************************
   ! check_btm2_top1 - spliting when TWO vertices on the btm
   !****************************************************************
   logical function check_btm2_top1 (hx1,postnCutSurf) result (lmsg)
      implicit none
      type(stl_tri), intent(inout)   :: hx1 
      real*8, intent(in)  :: postnCutSurf(3)
      real*8 :: tmpMaximum
      lmsg = .false.
      !  C  A | P | B
      tmpMaximum = max ( hx1%tvtx%vertexC(2),hx1%tvtx%vertexA(2) )
      if      ( tmpMaximum      < hx1%tvtx%vertexB(2)   ) then
         if   ( postnCutSurf(2) < hx1%tvtx%vertexB(2)   .and. &
                tmpMaximum      < postnCutSurf(2)       ) lmsg = .true. 
      end if
      !  B  C  | P | A
      tmpMaximum = max ( hx1%tvtx%vertexB(2), hx1%tvtx%vertexC(2) )
      if      ( tmpMaximum      < hx1%tvtx%vertexA(2)   ) then
         if   ( postnCutSurf(2) < hx1%tvtx%vertexA(2)   .and. &
                tmpMaximum      < postnCutSurf(2)       ) lmsg = .true. 
      end if
      !  A  B  | P | C
      tmpMaximum = max ( hx1%tvtx%vertexA(2),hx1%tvtx%vertexB(2) )
      if      ( tmpMaximum      < hx1%tvtx%vertexC(2)   ) then
         if   ( postnCutSurf(2) < hx1%tvtx%vertexC(2)   .and. &
                tmpMaximum      < postnCutSurf(2)       ) lmsg = .true. 
      end if
      return
   end function check_btm2_top1
   ! 
   !****************************************************************
   !  check_vertex_side_sign
   !****************************************************************
   logical function check_vertex_side_sign ( nmlvect, crossPt, posvect, sign_LR) result (msg) !#  
      ! If the output sign, {sign_LR}, is positive, it means that {posvect} is on the right or upper side
      ! of the plane passing through {crossPt}, having a normal directional vector, {nmlvect}. 
      ! For example, if the normal vector is to the right (+x) and crossPt is at the origin, then
      ! sign_LR would be positive if the posvector's x-position is positive. 
      ! use math 
      implicit none
      real*8, intent(in)   :: nmlvect (3), crossPt (3), posvect (3)
      real*8, intent(out)  :: sign_LR
      sign_LR      =  nmlvect .dot. (posvect - crossPt)
      if ( sign_LR > 0 ) then      ; msg = .true.
      else                         ; msg = .false.
      end if
      return
   end function check_vertex_side_sign   
   !*****************************************************************************************
   ! check_vertex_mat_side_signs ( nmvecCutSurf, postnCutSurf, vertexMat, sign_vec, iamorig)
   !*****************************************************************************************
   logical function check_vertex_mat_side_signs ( nmlvect,crossPt,pos_mat,sign_vec,iamorig) result (msg) !#  
      ! use math 
      implicit none
      real*8,  intent(in)   :: nmlvect (3), crossPt (3), pos_mat (3,3)
      real*8,  intent(inout)  :: sign_vec(3)
      integer, intent(out)  :: iamorig
      real*8  :: posvect(3)
      logical :: mymsg(3)
      msg       = .true.
      posvect = pos_mat(1,:) ; mymsg(1) = check_vertex_side_sign ( nmlvect, crossPt, posvect, sign_vec(1))
      posvect = pos_mat(2,:) ; mymsg(2) = check_vertex_side_sign ( nmlvect, crossPt, posvect, sign_vec(2))
      posvect = pos_mat(3,:) ; mymsg(3) = check_vertex_side_sign ( nmlvect, crossPt, posvect, sign_vec(3))
      iamorig =  whoisob (sign_vec) ;
      if (iamorig .eq. 0)      msg = .false.
      ! write(*,*) iamorig,msg 
      ! write(*,*)mymsg
      ! write(*,*)  sign_vec
   end function check_vertex_mat_side_signs
   !*******************************************************
   ! For an integer i, to calculate its self, plus1 and mnus1 in a cycle of Ndiv.
   ! 2018-07-29-09-02-52-AM
   ! 2018-07-29-10-40-47-AM 
   !*******************************************************
   !  imod_plus1 - returns a cyclic plus one
   !*******************************************************
   integer function imod_plus1 (i,Ndiv)  result (imodplus)
      implicit none 
      integer, intent(in) :: i, Ndiv
      imodplus = mod( Ndiv + i - 0 , Ndiv) + 1
   end function imod_plus1
   !***************************************************************
   ! imod_self0  - a dummy fuction that returns an integer itself. 
   !***************************************************************
   integer function imod_self0 (i,Ndiv)  result (imodself)
      implicit none 
      integer, intent(in) :: i, Ndiv
      imodself = mod( Ndiv + i - 1 , Ndiv) + 1
   end function imod_self0
   !*******************************************************
   ! imod_mnus1 - returns a cyclic minus one
   !*******************************************************
   integer function imod_mnus1 (i,Ndiv)  result (imodmnus)
      implicit none 
      integer, intent(in) :: i, Ndiv
      imodmnus = mod( Ndiv + i - 2 , Ndiv) + 1
   end function imod_mnus1
   !******************************************************************************
   ! isgn - an integer function that returns a sing of a real number: -1, 0, and 1. 
   !******************************************************************************
   integer function isgn (realval) result(nsign)
      implicit none
      real*8, intent(in)   :: realval
      character(128)       :: warning= &
           "'integer function isign (realval) result(nsign)' returns 0 for the sign of realval."
      if (realval > 0.0d0) then
         nsign = 1
      elseif (realval < 0.0d0) then
         nsign = -1
      elseif (realval == 0.0d0) then
         nsign = 0 ; 
         write(*,*) warning
      end if
   end function isgn
   !*******************************************************************
   !  whoisob - a wrapping function that select a vector element, 
   !            having a different sign from the two others.
   !*******************************************************************
   integer function whoisob (vector) result(iamob)
      implicit none
      real*8 vector(3)
      iamob = whichone (vector)
   end function whoisob
   !*******************************************************************
   !  whichone - the wrapped function of "whoisob"
   !*******************************************************************
   integer function whichone (vector) result(idiff)
      implicit none 
      integer :: i ! (needed when called in a mian ) imod_plus1,  imod_mnus1, imod_self0, isgn
      real*8  :: vector(3),     vecProd(3)
      integer :: isgn_self(3),  isgn_plus(3),  isgn_mnus(3),  isgn_prod(3)
      do i = 1, 3
         isgn_self(i) = isgn(vector(i) * 1.0d0                     )
         isgn_plus(i) = isgn(vector(i) * vector (imod_plus1 (i,3)) )         
         isgn_mnus(i) = isgn(vector(i) * vector (imod_mnus1 (i,3)) )
      enddo
      do i = 1, 3
         if ( (isgn_self(i) .eq. isgn_plus(i)) .and. (isgn_self(i) .eq. isgn_mnus(i)) ) idiff = i 
         if ( (isgn_self(i) .ne. isgn_plus(i)) .and. (isgn_self(i) .ne. isgn_mnus(i)) ) idiff = i 
      enddo
      if ( isgn_self(1) == isgn_self(2) .and.  & 
           isgn_self(2) == isgn_self(3) .and.  & 
           isgn_self(3) == isgn_self(1)) then
         ! write(*,*)  "All numbers have a same sign:", isgn_self(1)
         idiff = 0 
      end if
   end function whichone
   !*******************************************************************
   !  stl_tri_copy_matirxTOvertex
   !*******************************************************************
   logical function stl_tri_copy_matirxTOvertex (hx1) result(msg)
      implicit none 
      type(stl_tri), intent(inout)   :: hx1
      hx1%tvtx%vertexA          = hx1%tvtx%vertexMat(1,:)
      hx1%tvtx%vertexB          = hx1%tvtx%vertexMat(2,:)
      hx1%tvtx%vertexC          = hx1%tvtx%vertexMat(3,:)
      msg = .true.
   end function stl_tri_copy_matirxTOvertex
   !*******************************************************************
   ! 
   !*******************************************************************
   logical function stl_tri_copy_vertexTOmatirx (hx1) result(msg)
      implicit none 
      type(stl_tri), intent(inout)   :: hx1
      hx1%tvtx%vertexMat(1,:)   = hx1%tvtx%vertexA 
      hx1%tvtx%vertexMat(2,:)   = hx1%tvtx%vertexB 
      hx1%tvtx%vertexMat(3,:)   = hx1%tvtx%vertexC 
      msg = .true.
   end function stl_tri_copy_vertexTOmatirx
   !*******************************************************************
   ! 
   !*******************************************************************
   subroutine make_a_facet_association_list ( & 
        nout, nfacets, kfacet, hx1, RX, RY,RXref, RYref, Nrx, Nry, NRfacets,&
        Min_X,Max_X,Min_Y,Max_Y,gx,gy,LX,LY, assocList, msg, verboseOut )
      implicit none
      integer, intent(in)               :: nout, nfacets, kfacet
      type(stl_tri), intent(in)         :: hx1
      real*8 , intent(in)               :: RX, RY, RXref, RYref
      integer, intent(in)               :: Nrx, Nry, nrfacets
      integer, intent(inout)            :: assocList (Nrx,Nry,0:NRfacets)
      logical                           :: msg, verboseOut, myVerbose
      real*8            :: Min_X,Max_X,Min_Y,Max_Y,gx, gy, LX,LY
      real*8            :: FXmin, FYmin, FZmin
      real*8            :: FXmax, FYmax, FZmax
      integer           :: irgn_min, jrgn_min, irgn_max, jrgn_max
      integer           :: irgn, jrgn
      real*8            :: dmin1, dmax1
      real*8            :: vecX (3), vecY (3), vecZ (3)
      integer           :: mynout, nFacetCount
      integer           :: irun, jrun  
      external          :: min, max 
      ! vecX(1) = dmax1 (hx1%tvtx%vertexA(1) - RXref , gx/2.d0)
      ! vecX(2) = dmax1 (hx1%tvtx%vertexB(1) - RXref , gx/2.d0)
      ! vecX(3) = dmax1 (hx1%tvtx%vertexC(1) - RXref , gx/2.d0)
      ! vecX(1) = dmin1 (hx1%tvtx%vertexA(1) - RXref , LX)
      ! vecX(2) = dmin1 (hx1%tvtx%vertexB(1) - RXref , LX)
      ! vecX(3) = dmin1 (hx1%tvtx%vertexC(1) - RXref , LX)
      ! write(*,*) "Rx, Ry" ,  Rx, Ry
      ! Using the real values
      myVerbose = .false. 
      if(myVerbose) write(*,*)  "I am here."
      vecX(1) = hx1%tvtx%vertexA(1) ! - RXref 
      vecX(2) = hx1%tvtx%vertexB(1) ! - RXref 
      vecX(3) = hx1%tvtx%vertexC(1) ! - RXref       
      vecY(1) = hx1%tvtx%vertexA(2) ! - RYref
      vecY(2) = hx1%tvtx%vertexB(2) ! - RYref
      vecY(3) = hx1%tvtx%vertexC(2) ! - RYref
      if(myVerbose) write(*,*) vecX, vecY
      if(myVerbose) write(*,*)  "I am here again."
      ! making the min and max of real values 
      FXmin = dmin1 (vecX(1),vecX(2),vecX(3))  ; !FXmin = abs(FXmin)
      FXmax = dmax1 (vecX(1),vecX(2),vecX(3))  ; 
      FYmin = dmin1 (vecY(1),vecY(2),vecY(3))  ; !FYmin = abs(FYmin)
      FYmax = dmax1 (vecY(1),vecY(2),vecY(3))  ;
      if(myVerbose) write(*,*) FXmin, FXmax, FYmin, FYmax
      if(myVerbose) write(*,*)  "I am here again 3."
      irgn_min = INT (abs(FXmin - RXref ) / Rx ) + 1         ; 
      irgn_max = INT (   (FXmax - RXref ) / Rx ) + 1
      jrgn_min = INT (abs(FYmin - RYref ) / Ry ) + 1         ; 
      jrgn_max = INT (   (FYmax - RYref ) / Ry ) + 1
       ! write(123,*)irgn_min, irgn_max,jrgn_min, jrgn_max
      do irgn = irgn_min, irgn_max
         do jrgn = jrgn_min, jrgn_max
            ! nFacetCount = assocList(irgn, jrgn, 0)
            ! nFacetCount = nFacetCount + 1 
            nFacetCount                         = assocList(irgn, jrgn, 0) + 1
            assocList(irgn, jrgn, 0 )           = nFacetCount
            assocList(irgn, jrgn, nFacetCount ) = kfacet
            ! if ( nFacetCount .gt. NRfacets -1    ) then
            ! write(234,*)irgn, jrgn, irgn_min, irgn_max,jrgn_min, jrgn_max, NRfacets,&
            ! nFacetCount, assocList(irgn, jrgn, 0 ), kfacet;
            ! write(234,*)irgn, jrgn,  kfacet
                ! pause
            ! end if
         end do
      end do
      return
   end subroutine make_a_facet_association_list
   ! 
   !
   !       split_a_triangle_mat    : the same routine to split_a_triangle, but using matrix form.
   !       split_a_triangle        : the first version 
   !       split_a_triangle_v1     : an updated version using function check_vertex_side_sign
   !       check_vertex_side_sign
   !
   subroutine split_a_triangle_mat (ABCmat, crossPt, nmlvect, DEFmat, sgvec, ntri_tot, msg) !#  
      ! (colvec1, colvec2, colvec3, crossPt, nvector, Q12, Q23, Q31, sg1, sg2, sg3, ntri_tot, msg )
      ! 
      ! colvec1, colvec2, colvec3  : three vertex points of a triangle
      ! crossPt, nvector           : a point where a plane is made with normal vector, nvector
      ! message                    : if true, the the triangle is splitted, otherwise remain same. 
      ! ntri_tot                   : the number of total triangles after splitting, e.g.,
      !                                    1: no splitting, 4: regular splitting, and 
      !                                    2: degenerated splitting (any of two Q are equal.)
      ! Q12, Q23, Q31              : the three intersection positions with the plane along the three lines of the triangle.
      ! 2018-07-26-19-03-08-PM
                   !            C         !            A           !            B
                   !        F             !        D               !        E    
                   !    A   |   E         !    B   |   F           !    C   |   D
                   !        D             !        E               !        F    
                   !            B         !            C           !            A
                   !                      !                        !             
      ! use math
      implicit none
      real*8, intent(in)   :: ABCmat(3,3), crossPt(3), nmlvect(3)
      real*8, intent(out)  :: DEFmat(3,3)          ! 
      real*8, intent(out)  :: sgvec(3)             ! 
      logical,intent(out)  :: msg
      real*8               :: nvector(3)
      real*8               :: colvec1(3), colvec2(3), colvec3(3)
      real*8               :: Q12(3), Q23(3), Q31(3)
      real*8               :: sg1, sg2, sg3
      real*8               :: t12, t23, t31
      real*8               :: dvalnc, dvaln1, dvaln2, dvaln3
      real*8               :: testloc12, testloc23, testloc31
      integer              :: ntri_tot
      logical              :: lmsg
      ! inputs to call split_a_triangle
      colvec1 = ABCmat (1,:)  ;       
      colvec2 = ABCmat (2,:)  ;       
      colvec3 = ABCmat (3,:)  ;       
      nvector = nmlvect
      call split_a_triangle & 
           (colvec1, colvec2, colvec3, crossPt, nvector, Q12, Q23, Q31, sg1, sg2, sg3, ntri_tot, lmsg )
      ! outputs from  split_a_triangle
      DEFmat (1,:) = Q12 
      DEFmat (2,:) = Q23 
      DEFmat (3,:) = Q31 
      sgvec(1)     = sg1
      sgvec(2)     = sg2
      sgvec(3)     = sg3
      msg          = lmsg
   end subroutine split_a_triangle_mat
   ! 
   subroutine split_a_triangle & 
        (colvec1, colvec2, colvec3, crossPt, nvector, Q12, Q23, Q31, sg1, sg2, sg3, ntri_tot, msg )
      ! 
      ! colvec1, colvec2, colvec3  : three vertex points of a triangle
      ! crossPt, nvector           : a point where a plane is made with normal vector, nvector
      ! message                    : if true, the the triangle is splitted, otherwise remain same. 
      ! ntri_tot                   : the number of total triangles after splitting, e.g.,
      !                                    1: no splitting, 4: regular splitting, and 
      !                                    2: degenerated splitting (any of two Q are equal.)
      ! Q12, Q23, Q31              : the three intersection positions with the plane along the three lines of the triangle.
      ! 2018-07-26-19-03-08-PM
                   !            C         !            A           !            B
                   !        F             !        D               !        E    
                   !    A   |   E         !    B   |   F           !    C   |   D
                   !        D             !        E               !        F    
                   !            B         !            C           !            A
                   !                      !                        !             
      ! use math
      implicit none
      real*8, intent(in)   :: colvec1(3), colvec2(3), colvec3(3)
      real*8, intent(in)   :: crossPt(3), nvector(3)
      real*8, intent(out)  :: Q12(3), Q23(3), Q31(3)
      real*8, intent(out)  :: sg1, sg2, sg3
      logical,intent(out)  :: msg
      real*8               :: t12, t23, t31
      real*8               :: dvalnc, dvaln1, dvaln2, dvaln3
      real*8               :: testloc12, testloc23, testloc31
      integer              :: ntri_tot
      msg          = .false. 
      dvalnc       = nvector .dot. crossPt
      dvaln1       = nvector .dot. colvec1
      dvaln2       = nvector .dot. colvec2
      dvaln3       = nvector .dot. colvec3
      t12          =(dvalnc - dvaln1) / (dvaln1 - dvaln2)
      t23          =(dvalnc - dvaln2) / (dvaln2 - dvaln3)
      t31          =(dvalnc - dvaln3) / (dvaln3 - dvaln1)
      Q12          = colvec1 + t12 * (colvec1 - colvec2)
      Q23          = colvec2 + t23 * (colvec2 - colvec3)
      Q31          = colvec3 + t31 * (colvec3 - colvec1)
      testloc12    = (Q12 - colvec1) .dot. (Q12 - colvec2)
      testloc23    = (Q23 - colvec2) .dot. (Q23 - colvec3)
      testloc31    = (Q31 - colvec3) .dot. (Q31 - colvec1)
      if           ((testloc12 .lt. 0.0d0)   & 
           .and.   ( testloc23 .lt. 0.0d0) )  then
                                                Q31 = (colvec1 + colvec3) / 2.0d0 ; 
                                                msg = .true. ; ntri_tot=4
      elseif       ((testloc23 .lt. 0.0d0)   &
        .and.   ( testloc31 .lt. 0.0d0) ) then
                                                Q12 = (colvec1 + colvec2) / 2.0d0 ;
                                                msg = .true. ; ntri_tot=4 
      elseif       ((testloc31 .lt. 0.0d0)   &
        .and.  ( testloc12 .lt. 0.0d0) ) then
                                                Q23 = (colvec2 + colvec3) / 2.0d0 ; 
                                                msg = .true. ; ntri_tot=4
      end if
                                                sg1 = nvector .dot. (colvec1 - crossPt)
                                                sg2 = nvector .dot. (colvec2 - crossPt)
                                                sg3 = nvector .dot. (colvec3 - crossPt)
      if ( msg .eqv. .true. )  then
         if     ( magntd_vector(Q12 - Q23) == 0.0d0) then ; ntri_tot = 2 ;    
         elseif ( magntd_vector(Q23 - Q31) == 0.0d0) then ; ntri_tot = 2 ;   
         elseif ( magntd_vector(Q31 - Q12) == 0.0d0) then ; ntri_tot = 2 ;   
         end if
      end if
   end subroutine split_a_triangle
   ! 
   subroutine split_a_triangle_v1 & 
           (colvec1, colvec2, colvec3, crossPt, nvector, Q12, Q23, Q31, sg1, sg2, sg3, ntri_tot, msg )
      ! 
      ! colvec1, colvec2, colvec3  : three vertex points of a triangle
      ! crossPt, nvector           : a point where a plane is made with normal vector, nvector
      ! message                    : if true, the the triangle is splitted, otherwise remain same. 
      ! ntri_tot                   : the number of total triangles after splitting, e.g.,
      !                                    1: no splitting, 4: regular splitting, and 
      !                                    2: degenerated splitting (any of two Q are equal.)
      ! Q12, Q23, Q31              : the three intersection positions with the plane along the three lines of the triangle.
      ! 2018-07-26-19-03-08-PM
                      !            C         !            A           !            B
                      !        F             !        D               !        E    
                      !    A   |   E         !    B   |   F           !    C   |   D
                      !        D             !        E               !        F    
                      !            B         !            C           !            A
                      !                      !                        !             
      ! use mathstl
      ! use math
      implicit none
      real*8, intent(in)   :: colvec1(3), colvec2(3), colvec3(3)
      real*8, intent(in)   :: crossPt(3), nvector(3)
      real*8, intent(out)  :: Q12(3), Q23(3), Q31(3)
      real*8, intent(out)  :: sg1, sg2, sg3
      logical,intent(out)  :: msg
      real*8               :: t12, t23, t31
      real*8               :: dvalnc, dvaln1, dvaln2, dvaln3
      real*8               :: testloc12, testloc23, testloc31
      integer              :: ntri_tot
      ! logical              :: check_vertex_side_sign
      msg          = .false. 
      dvalnc       = nvector .dot. crossPt
      dvaln1       = nvector .dot. colvec1
      dvaln2       = nvector .dot. colvec2
      dvaln3       = nvector .dot. colvec3
      t12          =(dvalnc - dvaln1) / (dvaln1 - dvaln2)
      t23          =(dvalnc - dvaln2) / (dvaln2 - dvaln3)
      t31          =(dvalnc - dvaln3) / (dvaln3 - dvaln1)
      Q12          = colvec1 + t12 * (colvec1 - colvec2)
      Q23          = colvec2 + t23 * (colvec2 - colvec3)
      Q31          = colvec3 + t31 * (colvec3 - colvec1)
      testloc12    = (Q12 - colvec1) .dot. (Q12 - colvec2)
      testloc23    = (Q23 - colvec2) .dot. (Q23 - colvec3)
      testloc31    = (Q31 - colvec3) .dot. (Q31 - colvec1)
      if           ((testloc12 .lt. 0.0d0)   & 
           .and.   ( testloc23 .lt. 0.0d0) ) & 
           then
                                                   Q31 = (colvec1 + colvec3) / 2.0d0 ; 
                                                   msg = .true. ; ntri_tot=4
      elseif       ((testloc23 .lt. 0.0d0)   &
           .and.   ( testloc31 .lt. 0.0d0) ) &
           then
                                                   Q12 = (colvec1 + colvec2) / 2.0d0 ;
                                                   msg = .true. ; ntri_tot=4 
      elseif       ((testloc31 .lt. 0.0d0)   &
           .and.  ( testloc12 .lt. 0.0d0) )  & 
           then
                                                   Q23 = (colvec2 + colvec3) / 2.0d0 ; 
                                                   msg = .true. ; ntri_tot=4
      end if
      msg = check_vertex_side_sign ( nvector, crossPt, colvec1, sg1) ! sg1 = nvector .dot. (colvec1 - crossPt)
      msg = check_vertex_side_sign ( nvector, crossPt, colvec2, sg2) ! sg2 = nvector .dot. (colvec2 - crossPt)
      msg = check_vertex_side_sign ( nvector, crossPt, colvec3, sg3) ! sg3 = nvector .dot. (colvec3 - crossPt)
      if ( msg .eqv. .true. )  then
         if     ( magntd_vector(Q12 - Q23) == 0.0d0) then ; ntri_tot = 2 ;    
         elseif ( magntd_vector(Q23 - Q31) == 0.0d0) then ; ntri_tot = 2 ;   
         elseif ( magntd_vector(Q31 - Q12) == 0.0d0) then ; ntri_tot = 2 ;   
         end if
      end if
   end subroutine split_a_triangle_v1
   !


   !****************************************************************
   !  write_ascii_stl_all - to write an ascii stl file [2020-01-05-13-20-20-PM] 
   !****************************************************************
   subroutine write_ascii_stl_all (outputfile,hx,numtri,title)
      implicit none 
      ! integer,                  intent(in)      :: Nx,Ny
      character(len=4)     :: Cloop
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat, lin=10, i,j ,k
      integer :: ifct , numtri
      type(stl_tri),            intent(inout)   :: hx(numtri)
      character(len=128),       intent(inout)   :: outputfile, title
      numtri = size(hx)
      
      OPEN(lin, file=outputfile, status='replace')
         WRITE(lin,"('solid',2X,A128)") title 
         write(*,*) title
         ifct_loop: do ifct = 1, numtri 
            call write_ascii_stl_each (outputfile,lin,hx(ifct))
          enddo ifct_loop
         WRITE(lin,"('endsolid',2X,A128)") title 
      CLOSE(lin)
   end subroutine write_ascii_stl_all


   !****************************************************************
   !  write_ascii_stl_w_nnb_all - to write an ascii stl file 
   !****************************************************************
   subroutine write_ascii_stl_w_nnb_all (outputfile,hx,hxnnb,numtri,title)
      implicit none 
      ! integer,                  intent(in)      :: Nx,Ny
      character(len=4)     :: Cloop
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat, lin=10, i,j ,k
      integer :: ifct , numtri
      type(stl_tri),            intent(inout)   :: hx(numtri)
      type(stl_nnb),            intent(inout)   :: hxnnb(numtri)
      character(len=128),       intent(inout)   :: outputfile, title
      numtri = size(hx)
      
      OPEN(lin, file=outputfile, status='replace')
         WRITE(lin,"('solid',2X,A128)") title 
         write(*,*) title
         ifct_loop: do ifct = 1, numtri 
            call  write_ascii_stl_w_nnb_each (ifct,outputfile,lin,hx(ifct),hxnnb(ifct))
          enddo ifct_loop
         WRITE(lin,"('endsolid',2X,A128)") title 
      CLOSE(lin)
   end subroutine write_ascii_stl_w_nnb_all

   
   !****************************************************************
   !  write_ascii_stl_w_nnb_each - to write an ascii stl file [2020-01-05-13-21-43-PM] 
   !****************************************************************
   subroutine write_ascii_stl_w_nnb_each (ifct,outputfile,lin,hx1,hxnnb1)
      implicit none 
      type(stl_tri), intent(inout) :: hx1 
      type(stl_nnb), intent(inout) :: hxnnb1
      integer      , intent(in)    :: ifct
      character(len=128),       intent(inout)   :: outputfile
      integer   , intent(in) :: lin
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat,  i,j
      character(len=16)     :: mynumstring
      
      WRITE(lin,"(2x,'facet normal',1X,3(1x,ES16.8))")  hx1%tnvec
      WRITE(lin,"(4x,'outer loop')")
      WRITE(lin,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexA
      WRITE(lin,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexB
      WRITE(lin,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexC
      WRITE(lin,"(4x,'endloop')")
      WRITE(lin,"(4x,'endfacet',4x)",advance="no") 
      write(mynumstring,*) ifct
      WRITE(lin,"(A)",advance="no") trim(adjustl(mynumstring))
      write(mynumstring,*) hxnnb1%fidPair(1)    ; WRITE(lin,"(',',A)",advance="no") trim(adjustl(mynumstring))
      write(mynumstring,*) hxnnb1%fidPair(2)    ; WRITE(lin,"(',',A)",advance="no") trim(adjustl(mynumstring))
      write(mynumstring,*) hxnnb1%fidPair(3)    ; WRITE(lin,"(',',A)",advance="no") trim(adjustl(mynumstring))

      do j = 1, hxnnb1%vidPairMax
         write(mynumstring,*) hxnnb1%vidPair(j)  ;
         WRITE(lin,"(',',A)",advance="no") trim(adjustl(mynumstring))
      enddo
      WRITE(lin,*)
         
      
      ! WRITE(lin,"(4x,'endfacet',1X,30(1X,G10.0))") ifct, hxnnb1%fidPair
      ! , &
           ! (hxnnb1%vidPair(j),j=1,hxnnb1%vidPairMax)

      ! WRITE(lin,"(4x,'endfacet')") 
      ! WRITE(lin,"(4x,'endfacet',1X,30(1X,G10.0))") hxnnb1%fidPair, &
           ! (hxnnb1%vidPair(j),j=1,hxnnb1%vidPairMax)
      
      ! WRITE(lin,"(4x,'endfacet',1X,3(1X,I8))") hxnnb1%fidPair
   end subroutine write_ascii_stl_w_nnb_each

   subroutine write_ascii_stl_w_nnb_each_v0 (outputfile,lin,hx1,hxnnb1)
      implicit none 
      type(stl_tri), intent(inout) :: hx1 
      type(stl_nnb), intent(inout) :: hxnnb1
      
      character(len=128),       intent(inout)   :: outputfile
      integer   , intent(in) :: lin
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat,  i,j 
      WRITE(lin,"(2x,'facet normal',1X,3(1x,ES16.8))")  hx1%tnvec
      WRITE(lin,"(4x,'outer loop')")
      WRITE(lin,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexA
      WRITE(lin,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexB
      WRITE(lin,"(6x,'vertex',1X,3(1x,ES16.8))")  hx1%tvtx%vertexC
      WRITE(lin,"(4x,'endloop')")
      WRITE(lin,"(4x,'endfacet',1X,3(1X,G10.0))") hxnnb1%fidPair
      ! WRITE(lin,"(4x,'endfacet',1X,3(1X,I8))") hxnnb1%fidPair
   end subroutine write_ascii_stl_w_nnb_each_v0


   !****************************************************************
   ! cal_num_of_facets
   !****************************************************************
   integer function cal_num_of_facets (inputfile) result (numFacets)
      ! This function is to count the number of lines of inputfile.
      ! This works with only gfortran, and ifort give erroneous number. 
      character(200), intent(in)   :: inputfile
      character(200)               :: inputline
      integer                      :: eastat, lin, numSTL, mode7Facet
      lin = 10
      ! OPEN(lin, file=inputfile , status='old', action='read', position='rewind')
      OPEN(lin, file=inputfile , status='old')
      loop1: DO
         READ(lin,*,iostat=eastat) inputline
         IF (eastat < 0) THEN
            numvalues = numvalues + 1
            WRITE(*,*) trim(inputfile), ' :number of records =', numvalues-1
            numline   = numvalues - 1
            EXIT loop1
         ELSE IF (eastat > 0) THEN
            write(*,* ) "IOSTAT is ", eastat
            STOP 'IO-error'
         ENDIF
         numvalues = numvalues + 1
      END DO loop1
      CLOSE(lin)
      numSTL = numline - 2 
      mode7Facet = mod(numSTL,7)
      if (mode7Facet == 0 ) then
         numFacets   =  numSTL / 7
      else
         numFacets   = -1
      endif
      write(*,*)  " mod(numSTL,7) = ",   mode7Facet
      write(*,*)  " numFacets     = ",   numFacets
   end function cal_num_of_facets

! devel[2020-01-05-16-11-15-PM] 

   !****************************************************************
   !  read_ascii_stl_w_nnb_each - to write an ascii stl file [2020-01-05-13-21-43-PM] 
   !****************************************************************
   subroutine read_ascii_stl_w_nnb_each (entreefile,lin,hx1,hxnnb1)
      implicit none 
      type(stl_tri), intent(inout) :: hx1 
      type(stl_nnb), intent(inout) :: hxnnb1
      
      character(len=128),       intent(inout)   :: entreefile
      integer   , intent(in) :: lin
      character(len=4)     :: Cloop   
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat,  i,j 
      READ(lin,*) Cfacet, Cnormal, hx1%tnvec (1:3)
      READ(lin,*) Couter, Cloop   
      READ(lin,*) Cvertex,hx1%tvtx%vertexA (1:3)
      READ(lin,*) Cvertex,hx1%tvtx%vertexB (1:3)
      READ(lin,*) Cvertex,hx1%tvtx%vertexC (1:3)
      READ(lin,*) Cendloop
      READ(lin,*) Cendfacet, hxnnb1%fidPair, hxnnb1%vidPair(1:hxnnb1%vidPairMax)
   end subroutine read_ascii_stl_w_nnb_each
   
   !****************************************************************
   !  read_ascii_stl_w_nnb_all - to write an ascii stl file 
   !****************************************************************
   subroutine read_ascii_stl_w_nnb_all (entreefile,hx,hxnnb,numtri,title_stl)
      implicit none 
      ! integer,                  intent(in)      :: Nx,Ny
      character(len=4)     :: Cloop
      character(len=4)     :: Cfacet  
      character(len=5)     :: Couter  
      character(len=6)     :: Cnormal , Cvertex
      character(len=7)     :: Cendloop
      character(len=8)     :: Cendfacet
      integer :: numline, cal_num_of_file_lines ,eastat, lin=10, i,j ,k
      integer :: ifct , numtri
      type(stl_tri),            intent(inout)   :: hx(numtri)
      type(stl_nnb),            intent(inout)   :: hxnnb(numtri)
      character(len=128),       intent(inout)   :: entreefile, title_stl
      numtri = size(hx)
      
      OPEN(lin, file=entreefile, status='old')
         READ(lin,*) title_stl 
         ! write(*,*) title_stl
         ifct_loop: do ifct = 1, numtri
            call  read_ascii_stl_w_nnb_each (entreefile,lin,hx(ifct),hxnnb(ifct))
            hx(ifct)%tid = ifct
          enddo ifct_loop
         WRITE(lin,*)
      CLOSE(lin)
   end subroutine read_ascii_stl_w_nnb_all
   

   !****************************************************************
   !  read_ascii_nnbed_facet_all - to write an ascii stl file 
   !****************************************************************
   ! subroutine read_ascii_nnbed_facet_all (entreefile,nnbed_hx,numtri,title_stl)
   !    implicit none 
   !    ! integer,                  intent(in)      :: Nx,Ny
   !    character(len=4)     :: Cloop
   !    character(len=4)     :: Cfacet  
   !    character(len=5)     :: Couter  
   !    character(len=6)     :: Cnormal , Cvertex
   !    character(len=7)     :: Cendloop
   !    character(len=8)     :: Cendfacet
   !    integer :: numline, cal_num_of_file_lines ,eastat, lin=10, i,j ,k
   !    integer :: ifct , numtri
   !    type(stl_tri_nnb), intent(inout)   :: nnbed_hx(numtri)
   !    character(len=128),       intent(inout)   :: entreefile, title_stl
   !    numtri = size(hx)
      
   !    OPEN(lin, file=entreefile, status='old')
   !       READ(lin,*) title_stl 
   !       ! write(*,*) title_stl
   !       ifct_loop: do ifct = 1, numtri 
   !          call  read_ascii_stl_w_nnb_each (entreefile,lin,hx(ifct),hxnnb(ifct))
   !        enddo ifct_loop
   !       WRITE(lin,*)
   !    CLOSE(lin)
   ! end subroutine read_ascii_nnbed_facet_all
   
   

   
end module mathstl
! EOF 
