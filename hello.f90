PROGRAM Parallel_Hello_World
USE OMP_LIB


!$OMP PARALLEL

write(*,*)  'Hello from process: ', OMP_GET_THREAD_NUM()

!$OMP END PARALLEL

END
