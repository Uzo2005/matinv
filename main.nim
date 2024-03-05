import sequtils, strutils, algorithm
import fractionType, polynomialType
#TODOs
#implement eigen vectors and values
#orthogonal matrices
#diagonalisation of matrices

type 
    Matrix = object
        rows: int
        cols: int
        data: seq[Fraction]

    RowOperations = object
        intermediateStates: seq[Matrix]

    EigenPairs = object
        eigenValue: Fraction
        eigenVector: Matrix

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

proc `$`(mat: Matrix, tabIndents = 0): string =
    let indent = '\t'.repeat(tabIndents)
    result.add(indent & "[\n")

    for row in 0..<mat.rows:
        result.add(indent & "  [")
        for col in 0..<mat.cols:
            result.add(" " & $mat[row, col] & ",") 
        result.add("]\n")

    result.add(indent & "]\n")

proc toString(matrices: openArray[Matrix], tabIndents = 0): string =
    for index, matrix in matrices:
        result.add("(" & $(index+1) & ")\n" & `$`(matrix, tabIndents + 1) & "\n")

proc `$`(rowOperations: RowOperations): string =
    let 
        uniqueStates = rowOperations.intermediateStates.deduplicate
    result = $uniqueStates[^1]

proc showSteps(rowOperations: RowOperations): string =
    let 
        uniqueStates = rowOperations.intermediateStates.deduplicate
    result = uniqueStates.toString

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
    #this will set the elements of this matrix as data but will convert them to fractions

    result.rows = mat.rows
    result.cols = mat.cols
    result.data = data.mapit(it.toFraction)

proc withData(mat: Matrix, data: openArray[Fraction]): Matrix =
    #this will set the elements of this matrix as data

    result.rows = mat.rows
    result.cols = mat.cols
    result.data = data.toseq

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

proc determinant(mat: Matrix): Fraction {.raises: [IncompatibleDimensionsForMatrix].} =
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
    
proc adjoint(mat: Matrix): Matrix {.raises: [IncompatibleDimensionsForMatrix].} =
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

proc replaceColumn(mat: var Matrix, colIndex: int, data: openArray[Fraction]) =
    assert mat.rows == data.len

    for rowIndex, entry in data:
        mat[rowIndex, colIndex] = entry

proc replaceColumn(mat: var Matrix, colIndex: int, data: openArray[int]) =
    replaceColumn(mat, colIndex, data.mapit(it.toFraction))

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

proc findPivotColumn(rowData: seq[Fraction]): int =
    for index, data in rowData:
        if data != zeroFraction:
            return index
    
    return -1
    
proc makeZeroValueAtColumn(a, b: seq[Fraction], colIndex: int): seq[Fraction] =
    assert a.len == b.len, "How did your rows have different lengths?"

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
    # echo a, " ", b, " ", colIndex, " ", result

proc makeSurePivotsHaveTheirAbsoluteValue(mat: var Matrix) =
    for row in 0..<mat.rows:
        var rowData = mat.getRow(row)
        let pivotColumn = rowData.findPivotColumn

        if pivotColumn != -1:
            if rowData[pivotColumn] < 0:
                rowData.applyit(it * -1.toFraction)
                mat.replaceRow(row, rowData)

proc rowEchelonForm(mat: Matrix): RowOperations =
    #This will compute the row echelon form of a given matrix.
    #Remember that a matrix can have more than one echelon form

    var currentState = mat
    result.intermediateStates.add(currentState)
    var 
        pivotColumn: int
        pivotRow: int
    
    for rowIndex in 0..<mat.rows:
        let pivotColumnIndex = mat.getRow(rowIndex).findPivotColumn
        if pivotColumnIndex != -1:
            pivotRow = rowIndex
            pivotColumn = pivotColumnIndex
            break
    
    if pivotRow != 0:
        currentState.swapRows(pivotRow, 0)
        result.intermediateStates.add(currentState)
    
    var currentRowIndex: int
    while currentRowIndex < mat.rows:
        let 
            currentRowData = currentState.getRow(currentRowIndex)
            currentRowDataPivotColumn = currentRowData.findPivotColumn

        if currentRowDataPivotColumn == -1:
            if currentRowIndex != mat.rows - 1:
                var nextNonZeroRowIndex = -1
                for subsequentRowIndex in (currentRowIndex + 1)..<mat.rows:
                    let subsequentRowData = currentState.getRow(subsequentRowIndex)
                    if subsequentRowData.findPivotColumn != -1:
                        nextNonZeroRowIndex = subsequentRowIndex
                        break
                
                if nextNonZeroRowIndex != -1:
                    currentState.swapRows(currentRowIndex, nextNonZeroRowIndex)
                    continue
                else:
                    break
            else:
                break
        for subsequentRowIndex in (currentRowIndex + 1)..<mat.rows:
            let subsequentRowData = currentState.getRow(subsequentRowIndex)
            currentState.replaceRow(subsequentRowIndex, currentRowData.makeZeroValueAtColumn(subsequentRowData, currentRowDataPivotColumn))

        result.intermediateStates.add(currentState)
        inc currentRowIndex
    
    currentState.makeSurePivotsHaveTheirAbsoluteValue
    result.intermediateStates.add(currentState)
    result.intermediateStates = result.intermediateStates.deduplicate

proc reducedRowEchelonForm(mat: Matrix): RowOperations =
    #This will compute the reduced row echelon form of a given matrix

    var echelonForm = mat.rowEchelonForm.intermediateStates[^1]
    var currentRowIndex = echelonForm.rows - 1

    while currentRowIndex > 0:
        let 
            rowData = echelonForm.getRow(currentRowIndex)
            nonZeroElements = rowData.filterIt(it != zeroFraction)

        if nonZeroElements.len == 0:
            dec currentRowIndex
            continue

        let
            pivot = nonZeroElements[0]
            pivotIndex = rowData.find(pivot)
            scaledRowData = rowData.mapit(it / pivot)

        echelonForm.replaceRow(currentRowIndex, scaledRowData)

        for subsequentRowIndex in countDown(currentRowIndex - 1, 0):
            let subsequentRowData = echelonForm.getRow(subsequentRowIndex)
            echelonForm.replaceRow(subsequentRowIndex, scaledRowData.makeZeroValueAtColumn(subsequentRowData, pivotIndex))

        result.intermediateStates.add(echelonForm)
        dec currentRowIndex

    let 
        firstRow = echelonForm.getRow(0)
        firstRowNonZeroElements = firstRow.filterIt(it != zeroFraction)

    if firstRowNonZeroElements.len != 0:
        let
            firstRowPivot = firstRowNonZeroElements[0]
            scaledFirstRow = firstRow.mapit(it / firstRowPivot)
        echelonForm.replaceRow(0, scaledFirstRow)
        result.intermediateStates.add(echelonForm)
    
    echelonForm.data.applyit(it.toProperFraction)
    result.intermediateStates.add(echelonForm)
    result.intermediateStates = result.intermediateStates.deduplicate

    # result.data = echelonForm.data.mapit(it.toProperFraction)
    # result.rows = echelonForm.rows
    # result.cols = echelonForm.cols

proc solutionExists(augmentedMatrix: Matrix): bool =
    #Given the augmented Matrix of a system of linear equation, this will tell us if this system has a solution or not
    #Note that having a solution can mean that a unique solution exists or that the solutions are infinite

    let 
        rowEchelonForm = augmentedMatrix.rowEchelonForm.intermediateStates[^1]
        lastRowData = rowEchelonForm.getRow(rowEchelonForm.rows - 1)
        coefficientSumIsZero = lastRowData[0..^2].sum == zeroFraction
        
    result = true

    if coefficientSumIsZero and lastRowData[^1] != zeroFraction:
        return false

proc pivotColumn(row: openArray[Fraction]): int =
    #This will return the pivot column of a given row of a matrix (remember that the pivot column is the first column with a non-zero entry)
    #If no pivot column is found, -1 is returned
    for col in 0..<row.len:
        if row[col] != zeroFraction:
            return col

proc isInvertible(mat: Matrix): bool =
    #This checks if the given matrix is invertible
    assert mat.rows == mat.cols
    result = mat.reducedRowEchelonForm.intermediateStates[^1] == identityMat(mat.cols)

proc solveUsingRowReduction(augmentedMatrix: Matrix): string =
    #This takes the augmented matrix and outputs the solution if it exists
    #This will return "no solution exists" if no solution exists
    #This will segfault when the equation has greater than 26 unknowns, but I will cross that bridge when I reach it

    if not(solutionExists(augmentedMatrix)):
        return "No solution Exists For This System Of Linear Equations Because The Rightmost column Is A Pivot Column"

    let 
        reducedRowEchelonForm = augmentedMatrix.reducedRowEchelonForm
        allAlphabets = {'a'..'z'}.toseq
        neededVariables = allAlphabets[0..reducedRowEchelonForm.intermediateStates[^1].cols - 1]

    var basicVariables: seq[char]

    result.add("Solving for the values of " & join(neededVariables[0 .. ^2], ", ") & " and " & $neededVariables[^1] & " gives: \n")

    for rowIndex in countDown(augmentedMatrix.rows - 1, 0):
        var currentVariable: int
        let 
            rowData = reducedRowEchelonForm.intermediateStates[^1].getRow(rowIndex)

        basicVariables.add(neededVariables[rowData.pivotColumn])

        for column in 0..rowData.len - 2:
            if rowData[column] == zeroFraction:
                inc currentVariable
                continue
            
            let 
                coefficient = if rowData[column] < 0 : " - " & $(rowData[column] * - 1) elif(rowData[column] == 1.toFraction): "" else: " + " & $rowData[column]
            
            
            result.add(coefficient & neededVariables[currentVariable])
            inc currentVariable
        
        result.add(" = " & $rowData[^1] & "\n")
    
    result.add("The basic variables are " & basicVariables.reversed.join(", ") & " while the free variables are " & neededVariables.filterIt(it notin basicVariables).join(", "))

proc inverseFromRowOperations(mat: Matrix): tuple[rowOperations: RowOperations, inverse: Matrix] =
    #so to use row operations to find the inverse of the matrix,
    # I need to store all the operations that do the reduced row echelon form and then apply them to the Identity matrix in reverse
    #To do this, I just concat the identity Matrix with `mat`, do the row operations, and extract my inverse
    
    assert mat.cols == mat.rows
    result.inverse = initMat(mat.cols, mat.rows)

    result.rowOperations.intermediateStates.add(mat)

    let 
        identityMat = identityMat(mat.rows)

    var
        augmentedMatrix = initMat(mat.rows, mat.cols + mat.rows)
    
    for rowIndex in 0..<mat.rows:
        let augmentedRow = mat.getRow(rowIndex).concat(identityMat.getRow(rowIndex))
        augmentedMatrix.replaceRow(rowIndex, augmentedRow)
    
    result.rowOperations.intermediateStates.add(augmentedMatrix)
    
    let augmentedMatrixRowOperations = augmentedMatrix.reducedRowEchelonForm
    result.rowOperations.intermediateStates.add(augmentedMatrixRowOperations.intermediateStates)

    let augmentedMatrixContainingInverse = augmentedMatrixRowOperations.intermediateStates[^1]

    for rowIndex in 0..<mat.rows:
        let inverseRow = augmentedMatrixContainingInverse.getRow(rowIndex)[mat.cols .. ^1]
        result.inverse.replaceRow(rowIndex, inverseRow)
    result.rowOperations.intermediateStates.add(result.inverse)

proc LU(mat: Matrix): tuple[L: Matrix, U: Matrix] =
    #This will perform LU factorisation on `mat` and return the result as a `tuple[L, U]`
    var L = identityMat(mat.rows)
    var 
        rowEchelonForm = mat

    var 
        pivotColumn: int
        pivotRow: int

    for rowIndex in 0..<mat.rows:
        let pivotColumnIndex = mat.getRow(rowIndex).findPivotColumn
        if pivotColumnIndex != -1:
            pivotRow = rowIndex
            pivotColumn = pivotColumnIndex
            break

    if pivotRow != 0:
        rowEchelonForm.swapRows(pivotRow, 0)
    
    var currentRowIndex: int
    while currentRowIndex < mat.rows:
        let 
            currentRowData = rowEchelonForm.getRow(currentRowIndex)
        
        var
            currentRowNonZeroElements = currentRowData.filterIt(it != zeroFraction)
            currentRowPivotColumnIndex = -1

        if currentRowNonZeroElements.len > 0:
            currentRowPivotColumnIndex = currentRowData.find(currentRowNonZeroElements[0]) #TODO handle the case where pivot column is the last column

        if currentRowPivotColumnIndex != -1:
            let 
                pivotColumnData = newSeqWith(currentRowIndex, zeroFraction).concat rowEchelonForm.getColumn(currentRowPivotColumnIndex)[currentRowIndex .. ^1]
                pivotColumnDataForL = pivotColumnData.mapit(it / pivotColumnData[currentRowIndex])

            L.replaceColumn(currentRowIndex, pivotColumnDataForL)
        else:
            var nextNonZeroRowIndex = -1
            for subsequentRowIndex in (currentRowIndex + 1)..<mat.rows:
                let subsequentRowData = rowEchelonForm.getRow(subsequentRowIndex)
                if subsequentRowData.findPivotColumn != -1:
                    nextNonZeroRowIndex = subsequentRowIndex
                    break
            
            if nextNonZeroRowIndex != -1:
                rowEchelonForm.swapRows(currentRowIndex, nextNonZeroRowIndex)
                continue
            else:
                break

        
        for subsequentRowIndex in (currentRowIndex + 1)..<mat.rows:
            let subsequentRowData = rowEchelonForm.getRow(subsequentRowIndex)
            rowEchelonForm.replaceRow(subsequentRowIndex, currentRowData.makeZeroValueAtColumn(subsequentRowData, currentRowPivotColumnIndex))
        
        inc currentRowIndex 

    rowEchelonForm.makeSurePivotsHaveTheirAbsoluteValue

    result.L = L
    result.U = rowEchelonForm

proc characteristicEquation(mat: Matrix): PolynomialEquation =
    discard

proc useBackPropagationToObtainApproximateSolution(equation: PolynomialEquation): seq[Fraction] =
    discard

proc eigenPairs(mat: Matrix): seq[EigenPairs] =
    #This computes the eigenvalues and corresponding eigenvectors of matrix `mat` and returns all the eigenPairs as members of a sequence
    # If Ax = kx where k is a scalar called the eigenValue and x is a columnVector called the eigenVector, then (A - kI)x = 0. 
    # The characteristic equation of this matrix can be obtained by row reductions on the augmented matrix [(A - kI) 0]
    # If I can obtain the characteristic equation, I can obtain the approximate answers by backpropagation
    discard



when isMainModule:
    template generateData(row, cols: int): seq[int] = (1..row*cols).toseq
    
    # var mat2 = initMat(3, 3).withData([initFrac(13, 8), initFrac(-12, 8), initFrac(-1, 8), initFrac(-15, 8), initFrac(12, 8), initFrac(3, 8), initFrac(5, 4), zeroFraction, initFrac(-1, 4)])
    # var mat1 = initMat(3, 4).withData([0, 1, -4, 8, 2, -3, 2, 1, 5, -8, 7, 1])
    # var mat1 = initMat(3, 4).withData([1, -2, 1, 0, 0, 2, -8, 8, -4, 5, 9, -9])
    var mat1 = initMat(3, 3).withData([1, 1, 1, 0, 2, 3, 5, 5, 1])
    # var mat1 = initMat(4, 4).withData([3, 6, 4, 2, 8, 5, 3, 1, 7, 2, 3, 7, 8, 3, 2, 5])
    # var mat1 = initMat(3, 6).withData([0, 3, -6, 6, 4, -5, 3, -7, 8, -5, 8, 9, 3, -9, 12, -9, 6, 15])
    # var mat1 = initMat(3, 6).withData([3, -9, 12, -9, 6, 15, 0, 2, -4, 4, 2, -6, 0, 0, 0, 0, 1, 4])
    # var mat1 = initMat(3, 6).withData([0, 3, -6, 6, 4, -5, 3, -7, 8, -5, 8, 9, 3, -9, 12, -9, 6, 15])
    # var mat1 = initMat(3, 6).withData([1, 6, 2, -5, -2, -4, 0, 0, 2, -8, -1, 3, 0, 0, 0, 0, 1, 7])
    # var mat1 = initMat(3, 3).withData([1, 0, 0, 0, 1, 0, -4, 0, 1])
    # var mat1 = initMat(3, 3).withData([0, 1, 2, 1, 0, 3, 4, -3, 8])
    # var mat1 = initMat(5, 4).withData([2, -4, -2, 3, 6, -9, -5, 8, 2, -7, -3, 9, 4, -2, -2, -1, -6, 3, 3, 4])
    # var mat1 = initMat(4, 5).withData([2, 4, -1, 5, -2, -4, -5, 3, -8, 1, 2, -5, -4, 1, 8, -6, 0, 7, -3, 1])
    # var mat2 = initMat(5, 4).withData([2, -4, -2, 3, 6, -9, -5, 8, 2, -7, -3, 9, 4, -2, -2, -1, -6, 3, 3, 4])
    try: 
        # echo mat1.rowEchelonForm.showSteps
        # echo mat1.determinant
        # echo mat1.inverse
        # echo mat1.determinant
        echo mat1.inverseFromRowOperations.rowOperations.showSteps
        # echo mat1.inverse
        # echo mat1.reducedRowEchelonForm
        # echo mat1.solveUsingRowReduction
        # echo mat1.inverse
        # echo mat1.inverse * mat1
        # echo mat1.solutionExists
        # mat1.swapRows(1, 2)
        # echo mat1

        # let
        #     g = initFrac(3, 5)
        #     h = initFrac(1, 2)

        # echo lcm(g, h)

    except InverseOfSingularMatrix as e:
        echo e.msg
    
