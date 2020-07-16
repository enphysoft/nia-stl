module fufs
   ! frequently used format styles
   implicit none
   character(len=64) :: FMT_1RN         ="((2x,ES18.9))"
   character(len=64) :: FMT_4RN         ="(4(2x,ES18.9))"
   character(len=64) :: FMT_1C          ="(2X,A8,2X)"
   character(len=64) :: FMT_STRN        ="(A)"
   character(len=64) :: FMT_1C1I        ="(2X,A8,2X,I6)"
   character(len=64) :: FMT_1C4RN       ="(2X,A8,4(2x,ES18.9))"
   character(len=64) :: FMT_1S1I        ="((A8,1(2X, I8)))"
   character(len=64) :: FMT_1S8I        ="((A8,8(2X, I6)))"
   character(len=64) :: FMT_DTTM        ="((A8,8(I5)))"
   character(len=64) :: FMT_1S8I6       ="((A8,8(2X, I6)))"
   character(len=64) :: FMT_1S3S        ="((A8,3(2X,ES16.8)))"
   character(len=64) :: FMT_1S3E        ="((A8,3(2X, E16.8)))"
   character(len=64) :: FMT_1S3F        ="((A8,3(2X, F16.8)))"
   character(len=64) :: FMT_2I2RN       ="(2(2X,I5),2(2x,ES18.9))"
   character(len=64) :: FMT_3I3F        ="((3(2X,I6),3(2X, F16.8)))"
end module fufs
