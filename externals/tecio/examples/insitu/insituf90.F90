program insitu
    use iso_c_binding
    implicit none
#if defined TECIOMPI
    include "mpif.h"
#endif
    ! internal testing
    ! RUNFLAGS:none
    !
    ! Example FORTRAN program to write a partitioned binary in-situ
    ! SZPLT data file for Tecplot.
    !
    ! If TECIOMPI is #defined, this program may be executed with mpiexec with
    ! up to 3 MPI ranks (processes). In this case, it must be linked
    ! with the MPI version of TECIO.
    !

    include "tecio.f90"

    interface
        integer(c_int32_t) function initializeFile( &
#           if defined TECIOMPI
            mpiComm, &
#           endif
            numPartitions, fileHandle)
            use iso_c_binding
            implicit none
#           if defined TECIOMPI
            integer(c_int32_t), intent(in) :: mpiComm
#           endif
            integer(c_int32_t), intent(in) :: numPartitions
            type(c_ptr), intent(out) :: fileHandle
        end function initializeFile

        integer(c_int32_t) function createZone( &
            numPartitions, partitionOwners, fileHandle, zone)
            use iso_c_binding
            implicit none
            integer(c_int32_t), intent(in) :: numPartitions
            integer(c_int32_t), dimension(3), intent(in) :: partitionOwners
            type(c_ptr), intent(in) :: fileHandle
            integer(c_int32_t), intent(out) :: zone
        end function createZone

        integer(c_int32_t) function createData( &
#           if defined TECIOMPI
            commRank, &
#           endif
            numPartitions, partitionOwners, fileHandle, zone)
            use iso_c_binding
            implicit none
#           if defined TECIOMPI
            integer(c_int32_t), intent(in) :: commRank
#           endif
            integer(c_int32_t), intent(in) :: numPartitions
            integer(c_int32_t), dimension(3), intent(in) :: partitionOwners
            type(c_ptr), intent(in) :: fileHandle
            integer(c_int32_t), intent(in) :: zone
        end function createData

        integer(c_int32_t) function finalizeFile(fileHandle)
        use iso_c_binding
        type(c_ptr), intent(inout) :: fileHandle
        end function finalizeFile
    end interface

    integer(c_int32_t) :: returnValue = 0

    integer(c_int32_t), dimension(3) :: partitionOwners
    integer(c_int32_t) :: numPartitions = 3
    type(c_ptr) :: fileHandle = c_null_ptr
    integer(c_int32_t) :: zone

#   if defined TECIOMPI
    integer(c_int32_t) :: partition
    integer(c_int32_t) :: commSize
    integer(c_int32_t) :: commRank
    integer(c_int32_t) :: mainRank = 0
    integer(c_int32_t) :: ierr
    integer(c_int32_t) :: mpiComm = MPI_COMM_WORLD
    call MPI_Init(ierr)
    call MPI_Comm_size(mpiComm, commSize, ierr)
    call MPI_Comm_rank(mpiComm, commRank, ierr)
    do partition = 0, 2
        partitionOwners(partition + 1) = mod(partition, commSize)
    enddo
#   endif
    
    if (returnValue == 0) then
        returnValue = initializeFile( &
#               if defined TECIOMPI
            mpiComm, &
#               endif
            numPartitions, fileHandle)
    endif

    ! Create zone
    if (returnValue == 0) then
        returnValue = createZone(numPartitions, partitionOwners, fileHandle, zone)
    endif

    ! Create the connectivity and variable data
    if (returnValue == 0) then
        returnValue = createData( &
#               if defined TECIOMPI
            commRank, &
#               endif
            numPartitions, partitionOwners, fileHandle, zone)
    endif

    if (returnValue == 0) then
        returnValue = finalizeFile(fileHandle)
    endif
    
#   if defined TECIOMPI
    call MPI_Finalize(ierr)
#   endif

end program insitu

integer(c_int32_t) function initializeFile( &
#   if defined TECIOMPI
    mpiComm, &
#   endif
    numPartitions, fileHandle)
    use iso_c_binding
    implicit none
#   if defined TECIOMPI
    integer(c_int32_t), intent(in) :: mpiComm
#   endif
    integer(c_int32_t), intent(in) :: numPartitions
    type(c_ptr), intent(out) :: fileHandle
    
    integer(c_int32_t) :: fileFormat = 1
    integer(c_int32_t) :: debug = 1
    character(128) :: fileName
    character(24) :: varNames
    integer(c_int32_t) :: dataFileType = 0 ! FULL
    integer(c_int32_t) :: defaultVarType = 1 ! single precision
    type(c_ptr) :: gridFileHandle = c_null_ptr
    integer(c_int32_t) :: varSelection
    
#   if defined TECIOMPI
    integer(c_int32_t) :: mainRank = 0
#   endif

    include "tecio.f90"

    initializeFile = 0

    fileName = "insituf90.szplt"//char(0)
    varNames = "x y z p"//char(0)

    initializeFile = tecFileWriterOpen( &
        fileName, &
        "SIMPLE DATASET"//char(0), &
        varNames, &
        fileFormat, &
        dataFileType, &
        defaultVarType, &
        gridFileHandle, &
        fileHandle)
    
    if (initializeFile == 0) then
        initializeFile = tecFileSetDiagnosticsLevel(fileHandle, debug)
    endif
    
    if (initializeFile == 0) then
        ! Select a range of var 4 for output
        initializeFile = tecFileCreateVarSelection(fileHandle, 4, 180.0d0, 181.0d0, varSelection)
    endif

#   if defined TECIOMPI
    if (initializeFile == 0) then
        initializeFile = tecMPIInitialize(fileHandle, mpiComm, mainRank)
    endif
#   endif

    return
end

integer(c_int32_t) function createZone( &
    numPartitions, partitionOwners, fileHandle, zone)
    use iso_c_binding
    implicit none
    integer(c_int32_t), intent(in) :: numPartitions
    integer(c_int32_t), dimension(3), intent(in) :: partitionOwners
    type(c_ptr), intent(in) :: fileHandle
    integer(c_int32_t), intent(out) :: zone

    ! These should really be either declared globally or passed to the function
    integer(c_int32_t), parameter :: XDIM = 10
    integer(c_int32_t), parameter :: YDIM = 9
    integer(c_int32_t), parameter :: ZDIM = 8

    integer(c_int32_t) :: zoneType  = 5      ! Brick
    integer(c_int64_t) :: numNodes    = XDIM * YDIM * ZDIM ! Overall zone dimensions
    integer(c_int64_t) :: numCells    = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1)
    integer(c_int64_t) :: numFaceConnections   = 0
    integer(c_int32_t) :: shareConnectivityFromZone  = 0
    integer(c_int32_t) :: faceNeighborMode    = 0
    integer(c_int32_t), dimension(4) :: varTypes = (/1, 1, 1, 1/) ! 1 = real*4, 2 = real*8, 3 = integer*4, 4 = integer*2, 5 = CHARACTER
    integer(c_int32_t), dimension(4) :: shareVarFromZone = (/0, 0, 0, 0/) ! No variable sharing (have only one zone)
    integer(c_int32_t), dimension(4) :: valueLocations = (/1, 1, 1, 1/) ! 1 = Nodal, 0 = Cell-Centered
    integer(c_int32_t), dimension(4) :: passiveVarList = (/0, 0, 0, 0/) ! no passive variables
    character(1024) :: zoneTitle

    include "tecio.f90"

    createZone = 0

    zoneTitle = "partitioned Zone"//char(0)
    
    createZone = tecZoneCreateFE( &
        fileHandle, &
        zoneTitle, &
        zoneType, &
        numNodes, &
        numCells, &
        varTypes, &
        shareVarFromZone, &
        valueLocations, &
        passiveVarList, &
        shareConnectivityFromZone, &
        numFaceConnections, &
        faceNeighborMode, &
        zone)

#   if defined TECIOMPI
        ! Output partitions
        if (createZone == 0 .and. numPartitions .gt. 0) then
            createZone = tecZoneMapPartitionsToMPIRanks( &
                fileHandle, zone, numPartitions, partitionOwners)
        endif
#   endif

    return
end

integer(c_int32_t) function createData( &
#   if defined TECIOMPI
    commRank, &
#   endif
    numPartitions, partitionOwners, fileHandle, zone)
    use iso_c_binding
    implicit none
#   if defined TECIOMPI
    integer(c_int32_t), intent(in) :: commRank
#   endif
    integer(c_int32_t), intent(in) :: numPartitions
    integer(c_int32_t), dimension(3), intent(in) :: partitionOwners
    type(c_ptr), intent(in) :: fileHandle
    integer(c_int32_t), intent(in) :: zone

    ! These should really be either declared globally or passed to the function
    integer(c_int32_t), parameter :: XDIM = 10
    integer(c_int32_t), parameter :: YDIM = 9
    integer(c_int32_t), parameter :: ZDIM = 8

    integer(c_int64_t), dimension(3) :: iMin, iMax, jMin, jMax, kMin, kMax
    integer(c_int64_t), dimension(3) :: iDim, jDim, kDim
    integer(c_int64_t), dimension(3) :: pNNodes, pNCells ! Partition node and cell counts, including ghost items
    integer(c_int64_t) :: connectivityCount
    integer(c_int32_t) :: partition
    integer(c_int32_t) :: f
    integer(c_int32_t), dimension(:, :), allocatable :: connectivity
    real(c_float), dimension(:, :), allocatable ::  x, y, z, p
    integer :: allocateStatus
    integer :: maxCells, maxNodes
    integer :: i, j, k, index
    integer(c_int64_t), dimension(3) :: nGNodes, nGCells ! Partition ghost node and ghost cell counts
    integer(c_int32_t), dimension(:, :), allocatable, target :: ghostNodes, gNPartitions, gNPNodes, ghostCells
    ! Upper bound on counts of ghost nodes and cells
    integer :: maxGCells = XDIM * ZDIM + YDIM * ZDIM
    integer :: maxGNodes = 2 * (XDIM * ZDIM + YDIM * ZDIM)
    integer(c_int32_t) :: effectiveNumPartitions
    character(40) :: auxDataValue

    include "tecio.f90"

    interface
        subroutine GatherGhostNodesAndCells( &
            nGNodes, ghostNodes, gNPartitions, gNPNodes, &
            nGCells, ghostCells, iDim, jDim, kDim, jMin, jMax, maxGNodes, maxGCells)
            use iso_c_binding
            implicit none
            integer(c_int32_t), dimension(maxGNodes,3), target :: ghostNodes, gNPartitions, gNPNodes
            integer(c_int32_t), dimension(maxGCells,3), target :: ghostCells
            integer(c_int64_t), dimension(3) :: nGNodes, nGCells, iDim, jDim, kDim, jMin, jMax
            integer(c_int32_t), intent(in) :: maxGNodes, maxGCells
        end subroutine GatherGhostNodesAndCells
    end interface

    createData = 0

    if (numPartitions .gt. 0) then
        ! Divide the zone into 3 partitions, identified by the index ranges
        ! of an equivalent unpartitioned ordered zone.

        ! Partition 1 node range, which will include one layer of ghost cells on the IMAX boundary:
        iMin(1) = 0 ! We use zero-based indices here because it simplifies index arithmetic in C.
        iMax(1) = XDIM / 2 + 2
        jMin(1) = 0
        jMax(1) = YDIM
        kMin(1) = 0
        kMax(1) = ZDIM

        ! Partition 2; ghost cells on IMIN and JMAX boundaries:
        iMin(2) = iMax(1) - 3
        iMax(2) = XDIM
        jMin(2) = jMin(1)
        jMax(2) = YDIM / 2 + 2
        kMin(2) = kMin(1)
        kMax(2) = kMax(1)

        ! Partition 3; ghost cells on IMIN and JMIN boundaries:
        iMin(3) = iMin(2)
        iMax(3) = iMax(2)
        jMin(3) = jMax(2) - 3
        jMax(3) = YDIM
        kMin(3) = kMin(2)
        kMax(3) = kMax(2)

        effectiveNumPartitions = numPartitions
    else
        iMin(1) = 0
        iMax(1) = XDIM
        jMin(1) = 0
        jMax(1) = YDIM
        kMin(1) = 0
        kMax(1) = ZDIM

        iMin(2) = 0
        iMax(2) = 0
        jMin(2) = 0
        jMax(2) = 0
        kMin(2) = 0
        kMax(2) = 0

        iMin(3) = 0
        iMax(3) = 0
        jMin(3) = 0
        jMax(3) = 0
        kMin(3) = 0
        kMax(3) = 0

        effectiveNumPartitions = 1
    endif

    ! Local partition dimensions (of equivalent ordered zones)
    do partition = 1, 3
        if (partition .le. effectiveNumPartitions) then
            iDim(partition) = iMax(partition) - iMin(partition)
            jDim(partition) = jMax(partition) - jMin(partition)
            kDim(partition) = kMax(partition) - kMin(partition)
        else
            iDim(partition) = 0
            jDim(partition) = 0
            kDim(partition) = 0
        endif
    enddo

    ! Allocate memory for connectivity and variable values
    do partition = 1, effectiveNumPartitions
        if (partition .le. effectiveNumPartitions) then
            pNNodes(partition) = iDim(partition) * jDim(partition) * kDim(partition)
            pNCells(partition) = (iDim(partition) - 1) * (jDim(partition) - 1) * (kDim(partition) - 1)
        else
            pNNodes(partition) = 0
            pNCells(partition) = 0
        endif
    enddo
    
    maxCells = max(pNCells(1), pNCells(2), pNCells(3))
    maxNodes = max(pNNodes(1), pNNodes(2), pNNodes(3))
    allocate(connectivity(8 * maxCells, 3), stat = allocateStatus)
    if (allocateStatus == 0) allocate(x(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus == 0) allocate(y(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus == 0) allocate(z(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus == 0) allocate(p(maxNodes, 3), stat = allocateStatus)
    if (allocateStatus .ne. 0) stop 'Unable to allocate memory'
    
    ! Calculate variable and connectivity values for partitions.
    do partition = 1, effectiveNumPartitions
        ! Create variables
        do k = 0, kDim(partition) - 1
            do j = 0, jDim(partition) - 1
                do i = 0, iDim(partition) - 1
                    index = (k * jDim(partition) + j) * iDim(partition) + i + 1
                    x(index, partition) = real(i + iMin(partition) + 1)
                    y(index, partition) = real(j + jMin(partition) + 1)
                    z(index, partition) = real(k + kMin(partition) + 1)
                    p(index, partition) = real((i + iMin(partition) + 1) * (j + jMin(partition) + 1) * (k + kMin(partition) + 1))
                enddo
            enddo
        enddo

        ! connectivity
        do k = 0, kDim(partition) - 2
            do j = 0, jDim(partition) - 2
                do i = 0, iDim(partition) - 2
                    index = (k * (jDim(partition) - 1) + j) * (iDim(partition) - 1) + i

                    connectivity(8 * index + 1, partition) = (k * jDim(partition) + j) * iDim(partition) + i + 1
                    connectivity(8 * index + 2, partition) = &
                        connectivity(8 * index + 1, partition) + 1
                    connectivity(8 * index + 3, partition) = &
                        connectivity(8 * index + 1, partition) + iDim(partition) + 1
                    connectivity(8 * index + 4, partition) = &
                        connectivity(8 * index + 1, partition) + iDim(partition)
                    connectivity(8 * index + 5, partition) = &
                        connectivity(8 * index + 1, partition) + iDim(partition) * jDim(partition)
                    connectivity(8 * index + 6, partition) = &
                        connectivity(8 * index + 2, partition) + iDim(partition) * jDim(partition)
                    connectivity(8 * index + 7, partition) = &
                        connectivity(8 * index + 3, partition) + iDim(partition) * jDim(partition)
                    connectivity(8 * index + 8, partition) = &
                        connectivity(8 * index + 4, partition) + iDim(partition) * jDim(partition)
                enddo
            enddo
        enddo
    enddo

    if (numPartitions .gt. 0) then
        allocate(ghostNodes(maxGNodes, 3), stat = allocateStatus)
        if (allocateStatus == 0) allocate(gNPartitions(maxGNodes, 3), stat = allocateStatus)
        if (allocateStatus == 0) allocate(gNPNodes(maxGNodes, 3), stat = allocateStatus)
        if (allocateStatus == 0) allocate(ghostCells(maxGCells, 3), stat = allocateStatus)
        if (allocateStatus .ne. 0) stop 'Unable to allocate memory for ghost items'
        call GatherGhostNodesAndCells(nGNodes, ghostNodes, gNPartitions, gNPNodes, nGCells, ghostCells, &
            iDim, jDim, kDim, jMin, jMax, maxGNodes, maxGCells)
    endif

    do partition = 1, effectiveNumPartitions
#       if defined TECIOMPI
        if (numPartitions == 0 .or. partitionOwners(partition) == commRank) then
#       endif
        if (numPartitions .gt. 0) then
            if (createData == 0) &
                createData = tecFEPartitionCreate32( &
                    fileHandle, &
                    zone, &
                    partition, &
                    pNNodes(partition), &
                    pNCells(partition), &
                    nGNodes(partition), &
                    ghostNodes(1, partition), &
                    gNPartitions(1, partition), &
                    gNPNodes(1, partition), &
                    nGCells(partition), &
                    ghostCells(1, partition))
        endif

        if (createData == 0) then
            createData = tecZoneVarWriteFloatValues(fileHandle, zone, 1, partition, pNNodes(partition), x(1, partition))
        endif
        if (createData == 0) then
            createData = tecZoneVarWriteFloatValues(fileHandle, zone, 2, partition, pNNodes(partition), y(1, partition))
        endif
        if (createData == 0) then
            createData = tecZoneVarWriteFloatValues(fileHandle, zone, 3, partition, pNNodes(partition), z(1, partition))
        endif

        if (createData == 0) then
            createData = tecZoneVarWriteFloatValues(fileHandle, zone, 4, partition, pNNodes(partition), p(1, partition))
        endif

        ! Write out the connectivityCount
        connectivityCount = 8 * pNCells(partition)
        if (createData == 0) then
            createData = tecZoneNodeMapWrite32( &
                fileHandle, zone, partition, 1, connectivityCount, connectivity(1, partition))
        endif

#       if defined TECIOMPI
        endif
#       endif
    enddo

    if (numPartitions .gt. 0) then
        deallocate(ghostNodes, gNPartitions, gNPNodes, ghostCells)
    endif

    deallocate(x, y, z, p)

    return
end

integer(c_int32_t) function finalizeFile(fileHandle)
    use iso_c_binding
    include "tecio.f90"
    type(c_ptr), intent(inout) :: fileHandle
    finalizeFile = tecFileWriterClose(fileHandle)
    return
end

subroutine appendGhostItems( &
    nGhosts, ghosts, gPartitions, gPGhosts, partition, gIDim, gJDim, &
    gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd, oIDim, oJDim, &
    oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd, doPartition)
    use iso_c_binding
    implicit none
    integer(c_int64_t) :: nGhosts
    integer(c_int32_t) :: partition
    integer(c_int32_t), dimension(:) :: ghosts, gPartitions, gPGhosts
    integer(c_int64_t) :: gIDim, gJDim, gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd
    integer(c_int64_t) :: oIDim, oJDim, oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd
    logical :: doPartition

    integer :: i, j, k
    integer(c_int32_t) oI, oJ, oK
    do i = gIStart, gIEnd - 1
        do j = gJStart, gJEnd - 1
            do k = gKStart, gKEnd - 1
                nGhosts = nGhosts + 1
                ghosts(nGhosts) = (k * gJDim + j) * gIDim + i + 1
                if (doPartition) then
                    gPartitions(nGhosts) = partition
                    oI = i - gIStart + oIStart
                    oJ = j - gJStart + oJStart
                    oK = k - gKStart + oKStart
                    gPGhosts(nGhosts) = (oK * oJDim + oJ) * oIDim + oI + 1
                endif
            enddo
        enddo
    enddo
end


subroutine GatherGhostNodesAndCells( &
    nGNodes, ghostNodes, gNPartitions, gNPNodes, &
    nGCells, ghostCells, iDim, jDim, kDim, jMin, jMax, maxGNodes, maxGCells)
    use iso_c_binding
    implicit none
    integer(c_int32_t), target :: ghostNodes(maxGNodes,3), gNPartitions(maxGNodes, 3), gNPNodes(maxGNodes, 3), ghostCells(maxGCells, 3)
    integer(c_int64_t) :: nGNodes(3), nGCells(3), iDim(3), jDim(3), kDim(3), jMin(3), jMax(3)
    integer(c_int32_t), intent(in) :: maxGNodes, maxGCells
    !
    ! Assemble lists of ghost nodes and cells--nodes and cells near partition boundaries that
    ! coincide with those "owned" by neighboring partitions. For each set of coincident nodes
    ! or cells, exactly one partition must own the node or cell, and all other involved partitions
    ! must report it as a ghost node or cell.
    !
    ! Arbitrarily, we say that the first partition owns any nodes that do not overlay the interior
    ! of neighboring partitions. That is, it owns any nodes that its "real" (non-ghost) cells use.
    ! So only a single layer of nodes and cells--on its IMax boundary--are ghosts, and are owned
    ! by the second and third partitions. We use the same logic to assign ownership for nodes
    ! shared by partitions 2 and 3.
    !
    integer(c_int32_t), pointer, dimension(:) :: gnPtr, gnpPtr, gnpnPtr, gcPtr, unusedPtr
    interface
        subroutine appendGhostItems( &
        nGhosts, ghosts, gPartitions, gPGhosts, partition, gIDim, gJDim, &
        gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd, oIDim, oJDim, &
        oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd, doPartition)
        use iso_c_binding
        integer(c_int64_t) :: nGhosts
        integer(c_int32_t) :: partition
        integer(c_int32_t), dimension(:) :: ghosts, gPartitions, gPGhosts
        integer(c_int64_t) :: gIDim, gJDim, gIStart, gIEnd, gJStart, gJEnd, gKStart, gKEnd
        integer(c_int64_t) :: oIDim, oJDim, oIStart, oIEnd, oJStart, oJEnd, oKStart, oKEnd
        logical :: doPartition
        end subroutine appendGhostItems
    end interface

    nGNodes(1) = 0
    nGCells(1) = 0
    gnPtr => ghostNodes(1:maxGNodes, 1)
    gnpPtr => gnPartitions(1:maxGNodes, 1)
    gnpnPtr => gnpNodes(1:maxGNodes, 1)
    gcPtr => ghostCells(1:maxGCells, 1)
    ! Nodes owned by the second partition:
    call appendGhostItems( &
        nGNodes(1), gnPtr, gnpPtr, gnpnPtr, 2, &
        iDim(1), jDim(1), &                                 ! I- and J-dimensions
        iDim(1) - 1, iDim(1), 0, jDim(2) - 1, 0, kDim(1), & ! local index ranges
        iDim(2), jDim(2), &                                 ! I- and J-dimensions
        2, 3, 0, jDim(2) - 1, 0, kDim(2), .true.)           ! local index ranges
    ! Nodes owned by the third partition:
    call appendGhostItems( &
        nGNodes(1), gnPtr, gnpPtr, gnpnPtr, 3, &
        iDim(1), jDim(1), &                                       ! I- and J-dimensions
        iDim(1) - 1, iDim(1), jMin(3) + 2, jMax(3), 0, kDim(1), & ! local index ranges
        iDim(3), jDim(3), &                                       ! I- and J-dimensions
        2, 3, 2, jDim(3), 0, kDim(3), .true.)                     ! local index ranges
    ! Cells owned by the second partition:
    call appendGhostItems( &
        nGCells(1), gcptr, unusedPtr, unusedPtr, 2, &
        iDim(1) - 1, jDim(1) - 1, &                                 ! I- and J-dimensions
        iDim(1) - 2, iDim(1) - 1, 0, jDim(2) - 2, 0, kDim(1) - 1, & ! local index ranges
        iDim(2) - 1, jDim(2) - 1, &                                 ! I- and J-dimensions
        1, 2, 0, jDim(2) - 2, 0, kDim(2) - 1, .false.)              ! local index ranges
    ! Cells owned by the third partition:
    call appendGhostItems( &
        nGCells(1), gcptr, unusedPtr, unusedPtr, 3, &
        iDim(1) - 1, jDim(1) - 1, &                                           ! I- and J-dimensions
        iDim(1) - 2, iDim(1) - 1, jDim(2) - 2, jDim(1) - 1, 0, kDim(1) - 1, & ! local index ranges
        iDim(3) - 1, jDim(3) - 1, &                                           ! I- and J-dimensions
        1, 2, 1, jDim(3) - 1, 0, kDim(3), .false.)                            ! local index ranges

    nGNodes(2) = 0
    nGCells(2) = 0
    gnPtr => ghostNodes(1:maxGNodes, 2)
    gnpPtr => gnPartitions(1:maxGNodes, 2)
    gnpnPtr => gnpNodes(1:maxGNodes, 2)
    gcPtr => ghostCells(1:maxGCells, 2)
    ! Nodes owned by the first partition.
    call appendGhostItems( &
        nGNodes(2), gnPtr, gnpPtr, gnpnPtr, 1, &
        iDim(2), jDim(2), &                                       ! I- and J-dimensions
        0, 2, 0, jDim(2), 0, kDim(2), &                           ! local index ranges
        iDim(1), jDim(1), &                                       ! I- and J-dimensions
        iDim(1) - 3, iDim(1) - 1, 0, jDim(2), 0, kDim(1), .true.) ! local index ranges
    ! Nodes owned by the third partition.
    call appendGhostItems( &
        nGNodes(2), gnPtr, gnpPtr, gnpnPtr, 3, &
        iDim(2), jDim(2), &                             ! I- and J-dimensions
        2, iDim(2), jDim(2) - 1, jDim(2), 0, kDim(2), & ! local index ranges
        iDim(3), jDim(3), &                             ! I- and J-dimensions
        2, iDim(3), 2, 3, 0, kDim(3), .true.)           ! local index ranges
    ! Cells owned by the first partition.
    call appendGhostItems( &
        nGCells(2), gcptr, unusedPtr, unusedPtr, 1, &
        iDim(2) - 1, jDim(2) - 1, &                                        ! I- and J-dimensions
        0, 1, 0, jDim(2) - 1, 0, kDim(2) - 1, &                            ! local index ranges
        iDim(1) - 1, jDim(1) - 1, &                                        ! I- and J-dimensions
        iDim(1) - 3, iDim(1) - 2, 0, jDim(2) - 1, 0, kDim(1) - 1, .false.) ! local index ranges
    ! Cells owned by the third partition.
    call appendGhostItems( &
        nGCells(2), gcptr, unusedPtr, unusedPtr, 3, &
        iDim(2) - 1, jDim(2) - 1, &                                 ! I- and J-dimensions
        1, iDim(2) - 1, jDim(2) - 2, jDim(2) - 1, 0, kDim(2) - 1, & ! local index ranges
        iDim(3) - 1, jDim(3) - 1, &                                 ! I- and J-dimensions
        1, iDim(3) - 1, 1, 2, 0, kDim(3) - 1, .false.)              ! local index ranges

    nGNodes(3) = 0
    nGCells(3) = 0
    gnPtr => ghostNodes(1:maxGNodes, 3)
    gnpPtr => gnPartitions(1:maxGNodes, 3)
    gnpnPtr => gnpNodes(1:maxGNodes, 3)
    gcPtr => ghostCells(1:maxGCells, 3)
    ! Nodes owned by the first partition
    call appendGhostItems( &
        nGNodes(3), gnPtr, gnpPtr, gnpnPtr, 1, &
        iDim(3), jDim(3), &                                             ! I- and J-dimensions
        0, 2, 0, jDim(3), 0, kDim(3), &                                 ! local index ranges
        iDim(1), jDim(1), &                                             ! I- and J-dimensions
        iDim(1) - 3, iDim(1) - 1, jMin(3), jMax(3), 0, kDim(1), .true.) ! local index ranges
    ! Nodes owned by the second partition.
    call appendGhostItems( &
        nGNodes(3), gnPtr, gnpPtr, gnpnPtr, 2, &
        iDim(3), jDim(3), &                                       ! I- and J-dimensions
        2, iDim(3), 0, 2, 0, kDim(3), &                           ! local index ranges
        iDim(2), jDim(2), &                                       ! I- and J-dimensions
        2, iDim(2), jDim(2) - 3, jDim(2) - 1, 0, kDim(2), .true.) ! local index ranges
    ! Cells owned by the first partition
    call appendGhostItems( &
        nGCells(3), gcPtr, unusedPtr, unusedPtr, 1, &
        iDim(3) - 1, jDim(3) - 1, &                                              ! I- and J-dimensions
        0, 1, 0, jDim(3) - 1, 0, kDim(3) - 1, &                                  ! local index ranges
        iDim(1) - 1, jDim(1) - 1, &                                              ! I- and J-dimensions
        iDim(1) - 2, iDim(1) - 1, jMin(3), jMax(3) - 1, 0, kDim(1) - 1, .false.) ! local index ranges
    ! Cells owned by the second partition.
    call appendGhostItems( &
        nGCells(3), gcPtr, unusedPtr, unusedPtr, 2, &
        iDim(3) - 1, jDim(3) - 1, &                                        ! I- and J-dimensions
        1, iDim(3) - 1, 0, 1, 0, kDim(3) - 1, &                            ! local index ranges
        iDim(2) - 1, jDim(2) - 1, &                                        ! I- and J-dimensions
        1, iDim(2) - 1, jDim(2) - 2, jDim(2) - 1, 0, kDim(2) - 1, .false.) ! local index ranges
end
