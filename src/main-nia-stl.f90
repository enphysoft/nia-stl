! 
! Copyright © 2019-present, by Albert S. Kim
! Author: Albert S. Kim ( http://albertsk.org/ )
! Version: 1.0
! Package-Version: 20190630.0627
! Created: 06/30/2019
! Last modified: 07/15/2020
! Keywords: STL, nearest neighbor 
! Description: 
!  This program is, for a selected facet, to generate a list of 
! three edge-sharing, primary nearest facets, and
! multiple vertex-sharing, secondary nearest facets.   
! default entree = ascii-stl.stl
! output = ascii-stl_nia.stl, data/NNBfacet.dat, data/NNBindex.dat, data/NNBvertx.csv
!

program nia_stl
   use mathstl
   implicit none
   character (128) :: mystring
   character (128) :: entree_fileSTL
   character (128) :: defaultfileSTL="ascii-stl.stl"
   character (128) :: output_fileSTL="ascii-stl_nia.stl"
   character (128) :: output_fileNNB="data/NNBindex.dat"
   character (128) :: vertex_fileNNB="data/NNBvertx.csv"
   character (128) :: facets_fileNNB="data/NNBfacet.dat"
   character (128) :: title_stl
   character (128) :: prefix_nnbed="nia_", postfix_nnbed="_nia" 
   character (128) :: program_name
   character (128) :: args, arg1, arg1_dir, arg1_file, arg1_root, arg1_extn
   integer :: numtri, numline
   integer :: ipostn_sep, kpostn_sep,  ilenth_arg, ipostn_dot,  ilenth_file
   integer :: ix,jy, k
   integer :: ifct, kfct
   integer :: jvtx, lvtx
   integer :: Nfct, Nfct_o10, numPair_tmp
   integer :: fidPair_tmp(3), vtxPair_tmp(3), matPair_tmp(3,3)
   integer :: nargc
   real(8) :: vecA(3),vecB(3)
   logical :: sameVec
   integer      , allocatable :: numNborFacets(:), numNborVertices(:)
   type(stl_tri), allocatable :: myFacet(:), loadedFacet(:)
   type(stl_nnb), allocatable :: myNNbor(:), loadedNNbor(:)
   type(stl_tri_nnb) , allocatable :: nnbedFacet(:)
   
   nargc = iargc()
   call getarg(0,program_name)
   
   if ( nargc == 0 )  then
      entree_fileSTL=defaultfileSTL
   else if ( nargc == 1) then
      call getarg(1,arg1)
      entree_fileSTL=arg1
      ipostn_sep=scan(arg1,"/", back=.true.) 
      ilenth_arg= len(arg1)
      arg1_dir  = arg1(1:ipostn_sep)
      arg1_file = arg1(ipostn_sep+1:ilenth_arg)
      ipostn_dot=scan(arg1_file,".", back=.true.)
      ilenth_file= len(arg1_file)
      arg1_root = arg1_file(1:ipostn_dot-1)
      arg1_extn = arg1_file(ipostn_dot:ilenth_file)
      output_fileSTL=trim(adjustl(arg1_root))//trim(adjustl(postfix_nnbed))//trim(adjustL(arg1_extn))
      write(*,*) "Input  STL file: ", trim(adjustL(entree_fileSTL))
      write(*,*) "Output STL file: ", trim(adjustL(output_fileSTL))
   else
      write(*,*) " Correct usage: ", trim(adjustl(program_name))
      write(*,*) "            or: ", trim(adjustl(program_name)), " <stl-file-name>"
      stop
   end if
   
   write(*,"(' entree_fileSTL = ',A)") trim(entree_fileSTL)
   write(*,"(' output_fileSTL = ',A)") trim(output_fileSTL)
   write(*,"(' vertex_fileNNB = ',A)") trim(vertex_fileNNB)
   write(*,"(' facets_fileNNB = ',A)") trim(facets_fileNNB)
   write(*,"(' output_fileNNB = ',A)") trim(output_fileNNB)

   numline = cal_num_of_file_lines (entree_fileSTL)
   numtri  = (numline-2)/7
   Nfct    = numtri
   Nfct_o10= Nfct / 10
   
   write(*,*) "number of lines in input file = ", numline
   write(*,*) "number of triangles (facets)  = ", numtri
   
   allocate(myFacet(numtri), loadedFacet(numtri))
   allocate(myNNbor(numtri), loadedNNbor(numtri))
   allocate(nnbedFacet(numtri), numNborFacets(numtri), numNborVertices(numtri))
   
   call read_ascii_stl (entree_fileSTL,numtri,myFacet,title_stl)
   write(*,*) " reading STL file done. "
   
   do ifct = 1, Nfct
      myNNbor(ifct)%mvtx(1,1:3) = myFacet(ifct)%tvtx%vertexA(1:3)
      myNNbor(ifct)%mvtx(2,1:3) = myFacet(ifct)%tvtx%vertexB(1:3)
      myNNbor(ifct)%mvtx(3,1:3) = myFacet(ifct)%tvtx%vertexC(1:3)
   enddo
   
   numNborFacets   = 0
   numNborVertices = 0
   ifct_loop: do ifct = 1, Nfct
      kfct_loop: do kfct = ifct + 1 , Nfct
         numPair_tmp = 0
         matPair_tmp = 0
         fidPair_tmp = 0 
         vtxPair_tmp = 0             
         do jvtx = 1, 3
            do lvtx = 1, 3
               vecA = myNNbor(ifct)%mvtx(jvtx,1:3)
               vecB = myNNbor(kfct)%mvtx(lvtx,1:3)
               sameVec = chk_equal_vectors (vecA, vecB)
               if ( sameVec .eqv. .true. )  then
                  numPair_tmp = numPair_tmp + 1
                  fidPair_tmp (numPair_tmp) = kfct
                  matPair_tmp (jvtx, numPair_tmp) = lvtx
                  vtxPair_tmp (jvtx) = lvtx
               endif
            enddo
         enddo
         
         if (numPair_tmp == 1 ) then
             numNborVertices(ifct) = numNborVertices(ifct) + 1
             numNborVertices(kfct) = numNborVertices(kfct) + 1
             myNNbor(ifct)%VidPairMax = numNborVertices(ifct)
             myNNbor(kfct)%VidPairMax = numNborVertices(kfct)
             myNNbor(ifct)%VidPair(numNborVertices(ifct))  = kfct
             myNNbor(kfct)%VidPair(numNborVertices(kfct))  = ifct
          elseif (numPair_tmp == 2 ) then
             numNborFacets(ifct)                         = numNborFacets(ifct) + 1
             myNNbor(ifct)%fidPair(numNborFacets(ifct))  = kfct
             myNNbor(ifct)%matPair(numNborFacets(ifct),:)= vtxPair_tmp(:)
             numNborFacets(kfct)                         = numNborFacets(kfct) + 1
             myNNbor(kfct)%fidPair(numNborFacets(kfct))  = ifct
             myNNbor(kfct)%matPair(numNborFacets(kfct),:)= vtxPair_tmp(:)
         end if
         
      enddo kfct_loop

      if ( mod( ifct,Nfct_o10) == 0 )  then
         write(*,"(F6.2,'% done')") dble(ifct)/dble(Nfct) *100.d0
      else if (ifct == Nfct) then
         write(*,"(F6.2,'% done')") dble(ifct)/dble(Nfct) *100.d0
      end if
         
   enddo ifct_loop

   open(unit=21,file=vertex_fileNNB,status="replace")
   write(21,"('# ifct,fnnb1,fnnb2,fnnb3,vnnb1,vnnb2,vnnb...')")
   do ifct = 1, Nfct
      write(mystring,*) ifct                    ; write(21,"(    A)",advance="no") trim(adjustl(mystring))
      write(mystring,*) myNNbor(ifct)%FidPair(1); write(21,"(',',A)",advance="no") trim(adjustl(mystring))
      write(mystring,*) myNNbor(ifct)%FidPair(2); write(21,"(',',A)",advance="no") trim(adjustl(mystring))
      write(mystring,*) myNNbor(ifct)%FidPair(3); write(21,"(',',A)",advance="no") trim(adjustl(mystring))
      do k = 1, myNNbor(ifct)%VidPairMax
         write(mystring,*) myNNbor(ifct)%VidPair(k) ; write(21,"(',',A)",advance="no") trim(adjustl(mystring))
      enddo
      write(21,*)
   end do
   close(21)
   
   open(unit=22,file=output_fileNNB,status="replace")
   do ifct = 1, Nfct
      write(22,*) "neighbors of ", ifct, myNNbor(ifct)%fidPair(:)
      do ix = 1, 3
         write(22,*) "pairs of vtx ", ix, myNNbor(ifct)%matPair(:,ix)
      enddo
   enddo
   close(22)
   
   open(unit=23,file=facets_fileNNB,status="replace")
   write(23,"(A)") "solid nearest-neighbor-facets"
   do ifct = 1, Nfct
      write(23,"(A5,5X,4(2x,I8))") "facet", ifct, myNNbor(ifct)%fidPair(:)
      write(23,"(A4)") "pair" 
      do ix = 1, 3
         write(23,"(10X,4(2x,I8))") ix, myNNbor(ifct)%matPair(:,ix)
      enddo
      write(23,"(A7)") "endpair" 
      write(23,"(A8)") "endfacet" 
   enddo
   write(23,"(A)") "endsolid nearest-neighbor-facets"
   close(23)
   call write_ascii_stl_w_nnb_all (output_fileSTL,myFacet,myNNbor,numtri,title_stl)

   deallocate (myFacet, myNNbor, numNborFacets, loadedNNbor, loadedFacet, nnbedFacet )
   
 end program nia_stl
 
! EOF
