open (10, file='result_vectors/result_vector_32_test_1.data')
open (11, file='result_vectors/result_vector_64_test_1.data')

do i=1,64800
    read(10,*) a
    read(11,*) b
    if (abs(a-b) / abs(a).gt.1e-8) then
        write(6,*) i,a,b
    endif
enddo

end
