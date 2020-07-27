module m_HACApK_mine
	! a1(1:ndt, 1:kt)	a2(1:ndl, 1:kt)
	! nstrtt is leftCol of m
	! st_leafmtx is a2*a1 = Vm*Wm
	! st_leafmtx = st_leafmtxp%st_lf(ip) ?
	! zu is part of x
	! zbut = a1 * zu
	! zaut = a2 * zbut
	use m_HACApK_base
	implicit none
	
	interface length
		MODULE PROCEDURE lengthInt, lengthLogical
	endinterface length
	
	type :: newLeafmtxpParam
		integer*4 nlf ! number of sub-matrices
		integer*4 ktmax
		integer*4, pointer :: indexInOriginalLm(:)
		integer*4, pointer :: kt(:), ktOffset(:)
	endtype newLeafmtxpParam
	
contains
	integer function fullHcalcOrder(st_leafmtxp)
		type(st_HACApK_leafmtxp), intent(in) :: st_leafmtxp
		integer :: ip
		fullHcalcOrder = 0
		
		do ip = 1, st_leafmtxp%nlf
			fullHcalcOrder = fullHcalcOrder + leafCalcOrder(st_leafmtxp%st_lf(ip))
		enddo
	endfunction fullHcalcOrder
	
	
	integer function leafCalcOrder(st_leafmtx)
		type(sT_HACApK_leafmtx), intent(in) :: st_leafmtx
		if(st_leafmtx%ltmtx==1) then
			leafCalcOrder = (st_leafmtx%ndl + st_leafmtx%ndt) * st_leafmtx%kt
		elseif(st_leafmtx%ltmtx==2) then
			leafCalcOrder = st_leafmtx%ndl * st_leafmtx%ndt
		endif
	endfunction leafCalcOrder
	
	
	integer function lengthInt(array)
		integer, intent(in) :: array(:)
		integer, dimension(1) :: size
		size = shape(array)
		lengthInt = maxval(size)
	endfunction lengthInt
	
	
	integer function lengthLogical(array)
		logical, intent(in) :: array(:)
		integer, dimension(1) :: size
		size = shape(array)
		lengthLogical = maxval(size)
	endfunction lengthLogical
	
	
	logical function isAllInt(array, compare)
		integer, intent(in) :: array(:)
		integer, intent(in) :: compare
		integer :: i
		
		do i = 1, length(array)
			if(array(i) .ne. compare) then
				isAllInt = .FALSE.
				return
			endif
		enddo
		isAllInt = .TRUE.
	endfunction isAllInt
	
	
	logical function canSeparate(st_leafmtx, threads)
		type(st_HACApK_leafmtx), intent(in) :: st_leafmtx
		integer, intent(in) :: threads
		canSeparate = (st_leafmtx%kt >= threads .or. st_leafmtx%ndt >= threads)
	endfunction canSeparate
	
	
	logical function isAllTrue(array)
		logical, intent(in) :: array(:)
		integer :: i
		
		isAllTrue = .TRUE.
		do i = 1, length(array)
			isAllTrue = isAllTrue .and. array(i)
		enddo
	endfunction isAllTrue
	
	
	logical function isAllForAllProcs(array, compare, nProcs)
		include 'mpif.h'
		integer, intent(in) :: array(:)
		integer, intent(in) :: compare, nProcs
		integer :: ierr
		logical :: myIsAll
		logical, allocatable :: isAllArray(:)
		
		allocate(isAllArray(nProcs))
		myIsAll = isAllInt(array, compare)
		call MPI_AllGather(myIsAll, 1, MPI_LOGICAL, isAllArray, nProcs, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
		isAllForAllProcs = isAllTrue(isAllArray)
		
		deallocate(isAllArray)
	endfunction isAllForAllProcs
	
	
	logical function isAllProcs(element, compare, nProcs)
		include 'mpif.h'
		integer, intent(in) :: element, compare, nProcs
		integer :: ierr
		integer, allocatable :: elementArray(:)
!$omp barrier
		allocate(elementArray(nProcs))
		call MPI_AllGather(element, 1, MPI_INTEGER, elementArray, nProcs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
		print *, 'before isAllInt'
		isAllProcs = isAllInt(elementArray, compare)
		print *, 'after isAllInt'
		
		deallocate(elementArray)
	endfunction isAllProcs
	
	logical function isBigMatrix(st_leafmtx, HcalcOrder, threshold)
		type(st_HACApK_leafmtx), intent(in) :: st_leafmtx
		integer, intent(in) :: HcalcOrder
		real*8,  intent(in) :: threshold
		isBigMatrix = (leafCalcOrder(st_leafmtx) >= HcalcOrder * threshold)
	endfunction isBigMatrix
	
	logical function doesSpreadOverThreads(st_leafmtx, thWidth, threadsPerProc)
		type(st_HACApK_leafmtx), intent(in) :: st_leafmtx
		integer, intent(in) :: thWidth, threadsPerProc
		
		integer :: leftColID, rightColID
		leftColID = getThreadID(st_leafmtx%nstrtt, thWidth, threadsPerProc)
		rightColID = getThreadID(st_leafmtx%nstrtt + st_leafmtx%ndt - 1, thWidth, threadsPerProc)
		
		doesSpreadOverThreads = (st_leafmtx%ndt > thWidth .or. (leftColID .ne. rightColID))
	endfunction doesSpreadOverThreads
	
	integer function getThreadID(col, thWidth, threadsPerProc)
		integer, intent(in) :: col, thWidth, threadsPerProc
		getThreadID = mod(col / thWidth, threadsPerProc)
	endfunction getThreadID
	
	integer function getBroadestProcess(st_leafmtx, thWidth, threadsPerProc)
		type(st_HACApK_leafmtx), intent(in) :: st_leafmtx
		integer, intent(in) :: thWidth, threadsPerProc
		integer :: ndt, nstrtt, threadID, index, i
		integer, allocatable, dimension(:) :: vm
		
		ndt = st_leafmtx%ndt;  nstrtt = st_leafmtx%nstrtt
		allocate(vm(threadsPerProc))
		vm = 0
		
		do i = 1, ndt
			threadID = getThreadID(nstrtt + i - 1, thWidth, threadsPerProc)
			vm(threadID) = vm(threadID) + 1
		enddo
		
		index = getThreadID(nstrtt, thWidth, threadsPerProc)
		do while(vm(index) .ne. maxval(vm))
			index = mod(index + 1, threadsPerProc)
		enddo
		getBroadestProcess = index
		
		deallocate(vm)
	endfunction getBroadestProcess
	
	
	! thread starts from 0
	integer function loadBalanceQuot(size, threadsPerProc, thread)
		integer, intent(in) :: size, threadsPerProc
		integer, intent(in), optional ::  thread
		integer :: cantDivide
		if(mod(size, threadsPerProc) == 0) then
			loadBalanceQuot = size / threadsPerProc
		else
			cantDivide = mod(size, threadsPerProc)
			if(present(thread)) then
				if(thread < cantDivide) then
					loadBalanceQuot = size / threadsPerProc + 1
				else
					loadBalanceQuot = size / threadsPerProc
				endif
			else
				loadBalanceQuot = size / threadsPerProc + 1
			endif
		endif
	endfunction loadBalanceQuot
	
	
	integer function getRowReal8(matrix)
		real*8, intent(in) :: matrix(:,:)
		integer :: tempMatrix(2)
		
		tempMatrix = shape(matrix)
		getRowReal8 = tempMatrix(1)
	endfunction getRowReal8
	
	integer function getColReal8(matrix)
		real*8, intent(in) :: matrix(:,:)
		integer :: tempMatrix(2)
		
		tempMatrix = shape(matrix)
		getColReal8 = tempMatrix(2)
	endfunction getColReal8
	
	
	
	
	
	
	! method5
	! 小行列をつぶす
	subroutine makeNewLeafmtxp_test(st_leafmtxp, newLeafmtxp, div)
		type(st_HACApK_leafmtxp), intent(in)  :: st_leafmtxp
		type(st_HACApK_leafmtxp), intent(inout) :: newLeafmtxp
		integer, intent(in) :: div			!!	小行列を分割するときの幅
		
		type(st_HACApK_leafmtx) :: st_leafmtx
		type(st_HACApK_leafmtx), pointer :: destNLM
		integer :: nowNLMPindex, i, j, nlf, newNlf, newKt, ktStart, ktEnd, iExpand
		integer, allocatable, dimension(:) :: expandedSizeArray
		
		
		nowNLMPindex = 1
		nlf = st_leafmtxp%nlf
		allocate(expandedSizeArray(nlf))
		call makeNLMexpandedSizeArray(expandedSizeArray, st_leafmtxp, div)
		newNlf = sum(expandedSizeArray)
		
		newLeafmtxp%nd = st_leafmtxp%nd
		newLeafmtxp%nlf = newNlf
		newLeafmtxp%nlfkt = st_leafmtxp%nlfkt
		newLeafmtxp%ktmax = st_leafmtxp%ktmax
			if(newNlf > nlf)	newLeafmtxp%ktmax = div		!! 少しでも伸びていれば行列のktが分割されたことになるのでktmaxはdivになる
															!! divが10であれば、元のktmaxが9の時は新しいktmaxも9になるし、15だったら分割されて新しいktmaxは10になる　それ以上のktmaxは現れない
		allocate( newLeafmtxp%st_lf(newNlf) )
		
		
		if(.true.) then
			do i = 1, nlf
				st_leafmtx = st_leafmtxp%st_lf(i)
				
				if(st_leafmtx%ltmtx==2) then			! full matrix (not have to expand)
					newLeafmtxp%st_lf(nowNLMPindex) = st_leafmtx
				else
					do iExpand = 1, expandedSizeArray(i)
						destNLM => newLeafmtxp%st_lf( (nowNLMPindex + iExpand) - 1 )
						
						destNLM%ltmtx  = st_leafmtx%ltmtx
						destNLM%ndt    = st_leafmtx%ndt
						destNLM%ndl    = st_leafmtx%ndl
						
						destNLM%nstrtt = st_leafmtx%nstrtt
						destNLM%nstrtl = st_leafmtx%nstrtl
						
						if(iExpand < expandedSizeArray(i)) then		! ループの最後ではない場合
							newKt = div
							ktEnd = div * iExpand
						else										! ループの最後
							newKt = st_leafmtx%kt - div * (expandedSizeArray(i) - 1)
							ktEnd = st_leafmtx%kt
						endif
						
						destNLM%kt = newKt
						allocate( destNLM%a1(st_leafmtx%ndt, newKt), destNLM%a2(st_leafmtx%ndl, newKt) )
						
						ktStart = div * (iExpand - 1) + 1
						destNLM%a1(1:st_leafmtx%ndt, 1:newKt) = st_leafmtx%a1(1:st_leafmtx%ndt, ktStart:ktEnd)
						destNLM%a2(1:st_leafmtx%ndl, 1:newKt) = st_leafmtx%a2(1:st_leafmtx%ndl, ktStart:ktEnd)
						
						if(8 <= (nowNLMPindex + iExpand) - 1 .and. (nowNLMPindex + iExpand) - 1 <= 11) then
!							print *, "index", (nowNLMPindex + iExpand) - 1, "  ktStart", ktStart, "  ktEnd", ktEnd
						endif
					enddo
				endif
				nowNLMPindex = nowNLMPindex + expandedSizeArray(i)
			enddo
			
			if(.false.) then
				print *, "st_leafmtxp 6th a2    1st row   shape", shape(st_leafmtxp%st_lf(6)%a2)
				print *, st_leafmtxp%st_lf(6)%a2(1, :)
				print *, st_leafmtxp%st_lf(6)%a2(2, :)
				print *, st_leafmtxp%st_lf(6)%a2(3, :)
				
				print *, "newLeafmtxp 8,9,10,11th    1st row"
				print *, "shape 8th a2", shape(newLeafmtxp%st_lf(8)%a2)
				print *, newLeafmtxp%st_lf(8)%a2(1, :)
				print *, newLeafmtxp%st_lf(8)%a2(2, :)
				print *, newLeafmtxp%st_lf(8)%a2(3, :)
				print *, "shape 9th a2", shape(newLeafmtxp%st_lf(9)%a2)
				print *, newLeafmtxp%st_lf(9)%a2(1, :)
				print *, newLeafmtxp%st_lf(9)%a2(2, :)
				print *, newLeafmtxp%st_lf(9)%a2(3, :)
				print *, "shape 10th a2", shape(newLeafmtxp%st_lf(10)%a2)
				print *, newLeafmtxp%st_lf(10)%a2(1, :)
				print *, newLeafmtxp%st_lf(10)%a2(2, :)
				print *, newLeafmtxp%st_lf(10)%a2(3, :)
				print *, "shape 11th a2", shape(newLeafmtxp%st_lf(11)%a2)
				print *, newLeafmtxp%st_lf(11)%a2(1, :)
				print *, newLeafmtxp%st_lf(11)%a2(2, :)
				print *, newLeafmtxp%st_lf(11)%a2(3, :)
			endif
			
		endif
		
		deallocate(expandedSizeArray)
		print *, "end makeNewLeafmtxp_test"
	endsubroutine makeNewLeafmtxp_test
	
	
	
	
	
	
	! method6
	! つぶした後の小行列を計算し始める始点のみがわかるような構造体を作る
	subroutine makeNewLeafmtxpParam(st_leafmtxp, NLMparam, div)
		type(st_HACApK_leafmtxp), intent(in)  :: st_leafmtxp
		type(newLeafmtxpParam), intent(inout) :: NLMparam
		integer, intent(in) :: div			!!	小行列を分割するときの幅
		
		integer :: nowNLMPindex, index, i, j, nlf, newNlf
		integer, allocatable, dimension(:) :: expandedSizeArray, NLMindexArray
		integer, allocatable, dimension(:) :: sectionStart(:), sectionEnd(:), sectionKt(:)
		
		nlf = st_leafmtxp%nlf
		allocate(expandedSizeArray(nlf), NLMindexArray(nlf),   sectionStart(div), sectionEnd(div), sectionKt(div))
		
		call makeNLMexpandedSizeArray(expandedSizeArray, st_leafmtxp, div)
		call makeNLMexpandedLeafsIndexArray(expandedSizeArray, NLMindexArray, st_leafmtxp, div)
!		newNlf = NLMindexArray(nlf) + expandedSizeArray(nlf) - 1
		newNlf = sum(expandedSizeArray)
		
		NLMparam%nlf = newNlf
		NLMparam%ktmax = st_leafmtxp%ktmax
			if(newNlf > nlf)	NLMparam%ktmax = div		!! 少しでも伸びていれば行列のktが分割されたことになるのでktmaxはdivになる
		allocate( NLMparam%indexInOriginalLm(newNlf), NLMparam%kt(newNlf), NLMparam%ktOffset(newNlf) )
		
!		do i = 1, 50
		do i = 1, nlf
			nowNLMPindex = NLMindexArray(i)
			if( st_leafmtxp%st_lf(i)%ltmtx==2 ) then			! full matrix (not have to expand)
				NLMparam%kt(nowNLMPindex) = st_leafmtxp%st_lf(i)%kt
			else
				call setSectionStartEndKt(sectionStart, sectionEnd, sectionKt, st_leafmtxp%st_lf(i)%kt, div)
				
				if(omp_get_thread_num()==1)		print *, "i = ", i, "  expandedSizeArray = ", expandedSizeArray(i)
				
				do j = 1, expandedSizeArray(i)
					index = nowNLMPindex + j - 1
					NLMparam%indexInOriginalLm(index) = i
					NLMparam%kt(index) = sectionKt(j)
					NLMparam%ktOffset(index) = sectionStart(j)
					if(omp_get_thread_num()==1)		print '("originalIndex = ", i7, "  kt = ", i3, "  start = ", i3)',  NLMparam%indexInOriginalLm(index), NLMparam%kt(index), NLMparam%ktOffset(index)
				enddo
			endif
		enddo
		
!$omp master
		print '("newNlf = ", i7)', newNlf
!$omp end master
		
		deallocate(expandedSizeArray, NLMindexArray,   sectionStart, sectionEnd, sectionKt)
	endsubroutine makeNewLeafmtxpParam
	
	
	subroutine compareExpandedAndNLM(expandedSizeArray, NLMindexArray, max)
		integer, intent(in) :: expandedSizeArray(:), NLMindexArray(:), max
		integer :: i, counter
		counter = 0
		
		print *, "compareExpandedAndNLM"
		do i = 1, max
			if(expandedSizeArray(i) .ge. 2) then
				print '("i = ", i6, "  expandedSizeArray = ", i3, "  NLMindexArray = ", i6)', i, expandedSizeArray(i), NLMindexArray(i)
				counter = counter + 1
			endif
		enddo
	endsubroutine compareExpandedAndNLM
	
	
	! 各小行列を砕いたときそれぞれいくつに分割されるかの配列を作る
	! ltmtx==2のものは砕けない
	subroutine makeNLMexpandedSizeArray(expandedSizeArray, st_leafmtxp, div)
		integer, intent(out) :: expandedSizeArray(:)
		type(st_HACApK_leafmtxp), intent(in) :: st_leafmtxp
		integer, intent(in) :: div			!!	小行列を分割するときの幅
		
		integer :: i
		do i = 1, st_leafmtxp%nlf
			if(st_leafmtxp%st_lf(i)%ltmtx==2) then
				expandedSizeArray(i) = 1		! 分割できないものは1個として扱う
			else
				expandedSizeArray(i) = roundUp(st_leafmtxp%st_lf(i)%kt, div)
			endif
		enddo
	endsubroutine makeNLMexpandedSizeArray
	
	
	! 各小行列を砕いたものを新しい構造体に順に配置していったとき、
	! その構造体の中で元の小行列がそれぞれどのインデックスから始まるかを表す配列を作る
	subroutine makeNLMexpandedLeafsIndexArray(expandedSizeArray, NLMindexArray, st_leafmtxp, div)
		integer, intent(in) :: expandedSizeArray(:)
		integer, intent(out) :: NLMindexArray(:)
		type(st_HACApK_leafmtxp), intent(in) :: st_leafmtxp
		integer, intent(in) :: div			!!	小行列を分割するときの幅
		
		integer :: i, nowIndexInOldLF
		nowIndexInOldLF = 1
		
		do i = 1, st_leafmtxp%nlf
			NLMindexArray(i) = nowIndexInOldLF
			nowIndexInOldLF = nowIndexInOldLF + expandedSizeArray(i)
		enddo
	endsubroutine makeNLMexpandedLeafsIndexArray
	
	
	! newLeafmtxpは分割された低ランク行列で表されているが
	! newLeafmtxpの各インデックスの行列が元のleafmtxpの何番目のインデックスだったかを配列にする
	subroutine makeOriginalLfIndexArray(expandedSizeArray, NLM_originalIndexArray, originalNlf)
		integer, intent(in)  :: expandedSizeArray(:)
		integer, intent(out) :: NLM_originalIndexArray(:)
		integer, intent(in)  :: originalNlf			! st_leafmtxp%nlf
		
		integer :: i, j, indexInNLM
		indexInNLM = 1
		
		do i = 1, originalNlf
			do j = 1, expandedSizeArray(i)
				NLM_originalIndexArray(indexInNLM + j - 1) = i
			enddo
			indexInNLM = indexInNLM + expandedSizeArray(i)
		enddo
	endsubroutine makeOriginalLfIndexArray
	
	
	! 小行列を砕いて新しい構造体に格納する
	subroutine expandLeafIntoNLMp(st_leafmtx, newLeafmtxp, div, expandedSize, expandDestIndex)
		type(st_HACApK_leafmtx), intent(in) :: st_leafmtx
		type(st_HACApK_leafmtxp), intent(inout) :: newLeafmtxp
		integer, intent(in) :: div, expandedSize, expandDestIndex
		
		integer :: i, ndt, ndl, nowIndex
		integer, allocatable, dimension(:) :: sectionStart(:), sectionEnd(:), sectionKt(:)
		
		allocate(sectionStart(expandedSize), sectionEnd(expandedSize), sectionKt(expandedSize))
		call setSectionStartEndKt(sectionStart, sectionEnd, sectionKt, st_leafmtx%kt, div)
		
		do i = 1, expandedSize
			nowIndex = expandDestIndex + i - 1
			ndt = st_leafmtx%ndt
			ndl = st_leafmtx%ndl
			
			newLeafmtxp%st_lf(nowIndex)%ltmtx  = st_leafmtx%ltmtx
			newLeafmtxp%st_lf(nowIndex)%ndt    = ndt
			newLeafmtxp%st_lf(nowIndex)%ndl    = ndl
			newLeafmtxp%st_lf(nowIndex)%nstrtt = st_leafmtx%nstrtt
			newLeafmtxp%st_lf(nowIndex)%nstrtl = st_leafmtx%nstrtl
			newLeafmtxp%st_lf(nowIndex)%kt = sectionKt(i)
			allocate( newLeafmtxp%st_lf(nowIndex)%a1(ndt, sectionKt(i)) )
			allocate( newLeafmtxp%st_lf(nowIndex)%a2(ndl, sectionKt(i)) )
			
			newLeafmtxp%st_lf(nowIndex)%a1(1:ndt, :) = st_leafmtx%a1(1:ndt, sectionStart(i):sectionEnd(i))
			newLeafmtxp%st_lf(nowIndex)%a2(1:ndl, :) = st_leafmtx%a2(1:ndl, sectionStart(i):sectionEnd(i))
		enddo
		deallocate(sectionStart, sectionEnd, sectionKt)
	endsubroutine expandLeafIntoNLMp
	
	
	integer function roundUp(a, b)
		integer, intent(in) :: a, b
		if(mod(a, b)==0) then
			roundUp = a/b
		else
			roundUp = int(a/b) + 1
		endif
	endfunction roundUp
	
	
	! ある長さlengthをdivで配分していった時の開始点start、終点end、長さktをそれぞれ作る
	subroutine setSectionStartEndKt(start, end, kt, length, div)
		integer, intent(out) :: start(:), end(:), kt(:)
		integer, intent(in) :: length, div
		
		integer :: sections, sectionNo
		
		sections = roundUp(length, div)
		do sectionNo = 1, sections
			start(sectionNo) = (sectionNo - 1) * div + 1
			end(sectionNo)   = sectionNo * div				! 後で更新
			kt(sectionNo)    = div							! 後で更新
		enddo
		end(sections) = length								! end of last section
		kt(sections)  = length - start(sections) + 1		! kt of last section
	endsubroutine setSectionStartEndKt
	
	
	! method5で小行列を分割した後に計算量が均等になるように割り当てを行う
	subroutine methodFiveStaticOrderAllocation(newLeafmtxp, nThread, threadStartIndex, threadEndIndex)
		type(st_HACApK_leafmtxp), intent(in) :: newLeafmtxp
		integer, intent(in) :: nThread
		integer, intent(out) :: threadStartIndex(:), threadEndIndex(:)
		
		integer :: fullCalcOrder, threadCalcOrder, i, nowCalcOrder, nowIndex
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!! parallel region内で変数の宣言時に初期化するとその変数はthread privateにならない？
		!!! 上のfullCalcOrderについて、
		!!! integer :: fullCalcOrder と宣言するとthread private
		!!! integer :: fullCalcOrder = 0 と宣言するとthread sharedになる
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		fullCalcOrder = fullHcalcOrder(newLeafmtxp)
		threadCalcOrder = fullCalcOrder / nThread
		
		nowCalcOrder = 0
		nowIndex = 1
		
		! 先にthreadEndIndexを設定する
		do i = 1, newLeafmtxp%nlf
			nowCalcOrder = nowCalcOrder + leafCalcOrder(newLeafmtxp%st_lf(i))
			if(nowCalcOrder .ge. threadCalcOrder * nowIndex) then
				threadEndIndex(nowIndex) = i
				nowIndex = nowIndex + 1
			endif
			if(nowIndex == nThread) then
				threadEndIndex(nThread) = newLeafmtxp%nlf
				exit
			endif
		enddo
		
		threadStartIndex(1) = 1
		do i = 2, nThread
			threadStartIndex(i) = threadEndIndex(i-1) + 1
		enddo
	endsubroutine methodFiveStaticOrderAllocation
	
	! method9のreduction範囲を決める
	subroutine setLsLe(newLeafmtxp, nd, ls, le, startIndex, endIndex)
		type(st_HACApK_leafmtxp), intent(in) :: newLeafmtxp
		integer, intent(in) :: nd, startIndex, endIndex
		integer, intent(out) :: ls, le
		
		integer :: i, nstrtl, ndl
		ls = nd;		le = 1
		
		do i = startIndex, endIndex
			nstrtl = newLeafmtxp%st_lf(i)%nstrtl
			ndl = newLeafmtxp%st_lf(i)%ndl
			
			if(nstrtl < ls)				ls = nstrtl
			if(nstrtl + ndl - 1 > le)	le = nstrtl + ndl - 1
		enddo
		
!		print '("mythread  ", i2, "   start ", i7, "   end ", i7, "   ls ", i7, "   le ", i7)',   omp_get_thread_num(), startIndex, endIndex, ls, le
	endsubroutine setLsLe
endmodule m_HACApK_mine