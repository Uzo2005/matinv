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
    result.rows = mat.rows
    result.cols = mat.cols
    result.data = data.mapit(it.toFraction)

proc `*`(a, b: Matrix): Matrix =
    assert a.cols == b.rows, "Dimensions dont match"

    result = initMat(a.rows, b.cols)
    for row in 0..<a.rows:
        for col in 0..<b.cols:
            for value in 0..<a.cols:
                result[row, col] = result[row, col] + (a[row, value] * b[value, col])

proc `*`(a: Fraction, b: Matrix): Matrix =
    result = initMat(b.rows, b.cols)
    result.data = b.data.mapit(it*a)

proc cofactor(mat: Matrix, row, col: int): Matrix =
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
    result = initMat(mat.rows, mat.cols)

    for row in 0..<mat.rows:
        for col in 0..<mat.cols:
            result[col, row] = mat[row, col]

proc determinant(mat: Matrix): Fraction =
    if mat.shape == (2, 2):
        return (mat[0, 0] * mat[1, 1]) - (mat[0, 1] * mat[1, 0])
    
    result = initFrac()
    var currentMultiplier = initFrac(1)

    for col in 0..<mat.cols:
        let cofactorDet = determinant(mat.cofactor(0, col))
        result += currentMultiplier * mat[0, col] * cofactorDet
        currentMultiplier *= -1
    
proc adjoint(mat: Matrix): Matrix =
    #this does not check if this matrix even has a reasonable adjoint. 
    #this means it will compute adjoints for matrices which have a zero determinant

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

proc inverse(mat: Matrix): Matrix {.raises: [InverseOfSingularMatrix].} =
    if mat.shape == (2, 2):
        result = (determinant(mat).multiplicativeInverse) * adjoint(mat)

        return result

    let determinant = mat.determinant
    let adjoint = mat.adjoint

    if determinant == 0.0 :
        raise newException(InverseOfSingularMatrix, " This matrix does not have an inverse, because it has a determinant of Zero")

    result = (determinant.multiplicativeInverse) * adjoint


# when isMainModule:
#     template generateData(row, cols: int): seq[int] = (1..row*cols).toseq

#     let 
#         mat1 = initMat(4, 4).withData(@[1, -2, 2, 3, 2, -1, 6, 6, -2, 1, -4, -3, 1, -1, 4, 6])
#         mat2 = initMat(4, 1).withData(@[-7, -2, 0, -9])
#         invMat1 = inverse(mat1)

#     echo invMat1 * mat1
#     # # echo mat2
#     echo invMat1 * mat2




# when isMainModule:
#     template generateData(row, cols: int): seq[float] = (1..row*cols).toseq.mapit(it.float)
    
#     # let mat1 = initMat(4, 4).withData(@[3, 6, 4, 2, 8, 5, 3, 1, 7, 2, 3, 7, 8, 3, 2, 5].mapit(it.float))
#     # let mat1 = initMat(4, 4).withData(@[1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1].mapit(it.float))
#     # let mat1 = initMat(3, 3).withData(@[1, 1, 1, 0, 2, 3, 5, 5, 1].mapit(it.float))
#     # let mat1 = initMat(3, 3).withData(generateData(3, 3))
#     # let mat1 = initMat(4, 4).withData(@[1, 2, 2, -1, 2, -1, -1, 4, -1, 3, 4, -2, 3, 1, 1, -4].mapit(it.float))
#     # let mat1 = initMat(3, 3).withData(@[1, 0, 3, 1, 1, 0, 2, 0, 7].mapit(it.float))
#     # let mat1 = initMat(2, 2).withData(@[1, 2, 5, 3].mapit(it.float))
#     let mat1 = initMat(2, 2).withData(@[1, 2, 2, 3].mapit(it.float))
#     # let mat2 = initMat(4, 1).withData(@[5, 6, 5, 4].mapit(it.float))
#     let mat2 = initMat(2, 1).withData(@[3, 1].mapit(it.float))

#     try: 
#         echo mat1
#         # echo mat1.determinant
#         echo mat1.inverse
#         echo mat1.inverse * mat2

#         # echo mat2
#         # echo mat1.determinant
#         # echo mat1.adjoint
#         # echo (1/mat1.determinant) * (mat1 * mat1.adjoint)
#         # echo mat1.inverse * mat2
#         # echo mat1.inverse * mat1

#         # for data in (mat1.inverse * mat1).data:
#         #     echo data, " ", data.almostEqual(0.0)
#     except InverseOfSingularMatrix as e:
#         echo e.msg
    

