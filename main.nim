import sequtils
include fractionType
#TODOs
#implement cramer's rule
#implement eignen vectors

type Matrix = object
    rows: int
    cols: int
    data: seq[Fraction]

type InverseOfSingularMatrix = object of ValueError
type IncompatibleDimensionsForMatrix = object of ValueError

proc `[]`(mat: Matrix, row, col: int): Fraction =
    assert row < mat.rows, " Out of range value for row provided"
    assert col < mat.cols, " Out of range value for row provided"
    mat.data[(mat.cols * row) + col]

proc `[]=`(mat: var Matrix, row, col: int, val: Fraction) =
    assert row < mat.rows, " Out of range value for row provided"
    assert col < mat.cols, " Out of range value for row provided"
    mat.data[(mat.cols * row) + col] = val

proc `$`(mat: Matrix): string =
    result.add("[\n    ")

    for row in 0..<mat.rows:
        result.add("[")
        for col in 0..<mat.cols:
            result.add(" " & $mat[row, col] & ",") 
        result.add("]\n")

        if row != mat.rows - 1:
            result.add("    ")
        
    result.add("]\n")

proc initMat(rows, cols: int): Matrix =
    result.rows = rows
    result.cols = cols
    result.data = newSeqWith(rows * cols, initFrac())

proc identityMat(size: int) : Matrix =
    result = initMat(size, size)

    for index in 0..<size:
        result[index, index] = initFrac(1)

proc shape(mat: Matrix): tuple[rows: int, cols: int] =
    result.rows = mat.rows
    result.cols = mat.cols

proc withData(mat: Matrix, data: openArray[int]): Matrix =
    #this will set the elements of this matrix as data but will convert them to fraction

    result.rows = mat.rows
    result.cols = mat.cols
    result.data = data.mapit(it.toFraction)

proc `*`(a, b: Matrix): Matrix {.raises: IncompatibleDimensionsForMatrix.} =
    if a.cols != b.rows:
        raise newException(IncompatibleDimensionsForMatrix, "Matrix multiplication between A and B only works if A.cols == B.rows") 

    result = initMat(a.rows, b.cols)
    for row in 0..<a.rows:
        for col in 0..<b.cols:
            for value in 0..<a.cols:
                result[row, col] = result[row, col] + (a[row, value] * b[value, col])

proc `*`(a: Fraction, b: Matrix): Matrix =
    result = initMat(b.rows, b.cols)
    result.data = b.data.mapit(it*a)

proc cofactor(mat: Matrix, row, col: int): Matrix =
    #Given the position of an element in a matrix, this will return the cofactor of  the matrix at this position
    result = initMat(mat.rows - 1, mat.cols - 1)

    var currentIndex: int
    for rowIndex in 0..<mat.rows:
        if rowIndex == row:
            continue
        
        for colIndex in 0..<mat.cols:
            if colIndex == col:
                continue
            result.data[currentIndex] = mat[rowIndex, colIndex]
            inc currentIndex

proc transposed(mat: Matrix): Matrix =
    #This will return the transpose of a given matrix

    result = initMat(mat.rows, mat.cols)

    for row in 0..<mat.rows:
        for col in 0..<mat.cols:
            result[col, row] = mat[row, col]

proc determinant(mat: Matrix): Fraction {.raises: IncompatibleDimensionsForMatrix.} =
    #gives us the determinant of a given m x m Matrix

    if mat.rows != mat.cols:
        raise newException(IncompatibleDimensionsForMatrix, "Determinants only exist for square matrices")

    if mat.shape == (2, 2):
        return (mat[0, 0] * mat[1, 1]) - (mat[0, 1] * mat[1, 0])
    
    result = initFrac()
    var currentMultiplier = initFrac(1)

    for col in 0..<mat.cols:
        let cofactorDet = determinant(mat.cofactor(0, col))
        result += currentMultiplier * mat[0, col] * cofactorDet
        currentMultiplier *= -1
    
proc adjoint(mat: Matrix): Matrix {.raises: IncompatibleDimensionsForMatrix.} =
    #This will compute the adjoint of a given m x m matrix. By definition, the `adj(A) * A = xI` where x is a scalar

    #this does not check if this matrix even has a reasonable adjoint. 
    #this means it will compute adjoints for matrices which have a zero determinant!
    
    if mat.rows != mat.cols:
        raise newException(IncompatibleDimensionsForMatrix, "Adjoints only exist for square matrices")

    if mat.shape == (2, 2):
        result = initMat(2, 2)
        result[0, 0] = mat[1, 1]
        result[1, 1] = mat[0, 0]
        result[0, 1] = -1 * mat[0, 1]
        result[1, 0] = -1 * mat[1, 0]

        return result

    var adjointTransposed = initMat(mat.rows, mat.cols)

    for row in 0..<mat.rows:
        for col in 0..<mat.cols:
            let currentMultiplier = if (((row + col) and 1) == 0): initFrac(1) else: initFrac(-1)
            adjointTransposed[row, col] = currentMultiplier * determinant(mat.cofactor(row, col))

    result = adjointTransposed.transposed

proc inverse(mat: Matrix): Matrix {.raises: [InverseOfSingularMatrix, IncompatibleDimensionsForMatrix].} =
    #This will compute the inverse of a given matrix, it will raise `InverseOfSingularMatrix` exception when given a singular matrix

    if mat.shape == (2, 2):
        result = (determinant(mat).multiplicativeInverse) * adjoint(mat)

        return result

    let determinant = mat.determinant
    let adjoint = mat.adjoint

    if determinant == 0.0 :
        raise newException(InverseOfSingularMatrix, " This matrix does not have an inverse, because it has a determinant of Zero")

    result = (determinant.multiplicativeInverse) * adjoint

proc appendRow(mat: var Matrix, data: openArray[Fraction]) =
    #This will append `data` as the last row of `mat`
    assert mat.cols == data.len
    mat.data.add(data)
    mat.rows = mat.rows + 1

proc appendRow(mat: var Matrix, data: openArray[int]) =
    appendRow(mat, data.mapit(it.toFraction))

proc replaceRow(mat: var Matrix, rowIndex: int, data: openArray[Fraction]) =
    #This will replace the row at `rowIndex` of `mat` with `data`
    assert mat.cols == data.len
    mat.data = mat.data[0..(mat.cols * rowIndex ) - 1].concat(data.toseq).concat(mat.data[(mat.cols * rowIndex ) + mat.cols .. ^1])

proc replaceRow(mat: var Matrix, rowIndex: int, data: openArray[int]) =
    replaceRow(mat, rowIndex, data.mapit(it.toFraction))

proc getColumn(mat: Matrix, colIndex: int): seq[Fraction] =
    #This returns the `colIndex` column of `mat`

    for row in 0..<mat.rows:
        result.add(mat[row, colIndex])

proc getRow(mat: Matrix, rowIndex: int): seq[Fraction] =
    #This returns the `rowIndex` row of `mat`
    result = mat.data[(mat.cols * rowIndex) .. (mat.cols * rowIndex) + mat.cols - 1]

proc isNonZero(a: openArray[Fraction]): bool =
    #Given the column or row of a matrix, this will tell us if it has any non-zero element
    a.filterIt(it != zeroFraction).len > 0

proc swapRows(mat: var Matrix, a, b: int) =
    #This will swap rows `a` and `b` of `mat`

    let 
        rowA = mat.getRow(a)
        rowB = mat.getRow(b)
    mat.replaceRow(a, rowB)
    mat.replaceRow(b, rowA)

proc rowEchelonForm(mat: Matrix): Matrix =
    #This willl compute the row echelon form of a given matrix.
    #Remember that a matrix can have more than one echelon form

    proc makeZeroValueAtColumn(a, b: seq[Fraction], colIndex: int): seq[Fraction] =
        assert a.len == b.len, "How did your rows have different lengths"
        if a[colIndex] == zeroFraction or b[colIndex] == zeroFraction:
            return b
        let 
            lcmAtIndex = lcm(a[colIndex], b[colIndex])
            aMultiplier = lcmAtIndex / a[colIndex]
            bMultiplier = lcmAtIndex / b[colIndex]
            scaledA = a.mapit(it * aMultiplier)
            scaledB = b.mapit(it * bMultiplier)

        for index in 0..<a.len:
            result.add(scaledA[index] - scaledB[index])

    result = mat

    var 
        pivotColumn: int
        pivotRow: int

    for colIndex in 0..<mat.cols:
        if mat.getColumn(colIndex).isNonZero:
            pivotColumn = colIndex
            break

    for rowIndex in 0..<mat.rows:
        if mat[rowIndex, pivotColumn] != zeroFraction:
            pivotRow = rowIndex
            break

    if pivotRow != 0:
        result.swapRows(pivotRow, 0)
    
    var currentRowIndex: int
    while currentRowIndex < mat.rows:
        let currentRowData = result.getRow(currentRowIndex)


        for subsequentRowIndex in (currentRowIndex + 1)..<mat.rows:
            let subsequentRowData = result.getRow(subsequentRowIndex)
            result.replaceRow(subsequentRowIndex, currentRowData.makeZeroValueAtColumn(subsequentRowData, currentRowIndex))

        inc currentRowIndex
    

proc reducedRowEchelonForm(mat: Matrix): Matrix =
    #This will compute the reduced row echelon form of a given matrix
    proc makeZeroValueAtColumn(a, b: seq[Fraction], colIndex: int): seq[Fraction] =
        assert a.len == b.len, "How did your rows have different lengths"
        if a[colIndex] == zeroFraction or b[colIndex] == zeroFraction:
            return b
        let 
            lcmAtIndex = lcm(a[colIndex], b[colIndex])
            aMultiplier = lcmAtIndex / a[colIndex]
            bMultiplier = lcmAtIndex / b[colIndex]
            scaledA = a.mapit(it * aMultiplier)
            scaledB = b.mapit(it * bMultiplier)

        for index in 0..<a.len:
            result.add(scaledA[index] - scaledB[index])

    var echelonForm = mat.rowEchelonForm

    var currentRowIndex = echelonForm.rows - 1

    while currentRowIndex > 0:
        let 
            rowData = echelonForm.getRow(currentRowIndex)
            pivot = rowData.filterIt(it != zeroFraction)[0]
            pivotIndex = rowData.find(pivot)
            scaledRowData = rowData.mapit(it / pivot)

        echelonForm.replaceRow(currentRowIndex, scaledRowData)

        for thisRowIndex in countDown(currentRowIndex - 1, 0):
            let thisRowData = echelonForm.getRow(thisRowIndex)
            echelonForm.replaceRow(thisRowIndex, scaledRowData.makeZeroValueAtColumn(thisRowData, pivotIndex))

        dec currentRowIndex

    let 
        firstRow = echelonForm.getRow(0)
        firstRowPivot = firstRow.filterIt(it != zeroFraction)[0]
        scaledFirstRow = firstRow.mapit(it / firstRowPivot)
    
    echelonForm.replaceRow(0, scaledFirstRow)
    result.data = echelonForm.data.mapit(it.toProperFraction)
    result.rows = echelonForm.rows
    result.cols = echelonForm.cols

proc solutionExists(augmentedMatrix: Matrix): bool =
    #Given the augmented Matrix of a system of linear equation, this will tell us if this system has a solution or not
    #Note that having a solution can mean that a unique solution exists or that the solutions are infinite

    let reducedRowEchelonForm = augmentedMatrix.reducedRowEchelonForm

    for rowIndex in 0..<augmentedMatrix.rows:
        let 
            rowData = augmentedMatrix.getRow(rowIndex)
            coefficientSumIsZero = rowData[0..^2].sum == zeroFraction

        result = true

        if coefficientSumIsZero and rowData[^1] != zeroFraction:
            return false





when isMainModule:
    template generateData(row, cols: int): seq[int] = (1..row*cols).toseq
    
    # var mat1 = initMat(4, 4).withData([3, 6, 4, 2, 8, 5, 3, 1, 7, 2, 3, 7, 8, 3, 2, 5])
    # var mat1 = initMat(3, 6).withData([0, 3, -6, 6, 4, -5, 3, -7, 8, -5, 8, 9, 3, -9, 12, -9, 6, 15])
    var mat1 = initMat(3, 6).withData([1, 6, 2, -5, -2, -4, 0, 0, 2, -8, -1, 3, 0, 0, 0, 0, 1, 7])
   
    try: 
        echo mat1
        echo mat1.reducedRowEchelonForm
        echo mat1.solutionExists
        # mat1.swapRows(1, 2)
        # echo mat1

        # let
        #     g = initFrac(3, 5)
        #     h = initFrac(1, 2)

        # echo lcm(g, h)

    except InverseOfSingularMatrix as e:
        echo e.msg
    

