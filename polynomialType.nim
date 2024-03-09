import strformat
import ./fractionType

type 
    InfixOperations = enum
        Plus = "+"
        Minus = "-"
        Multiplication = "*"
        Division = "/"

    PolynomialTerm* = object
        case isMultipleTerms: bool
            of true:
                variables: seq[PolynomialTerm]
                infixOperations: seq[InfixOperations]
            else:
                coefficient: Fraction
                variable: char
        degree: Fraction
                # derivative: Fraction

    PolynomialEquation* = object #all polynomial equations implicitly equal 0
        variables: seq[PolynomialTerm] 
        infixOperations: seq[InfixOperations]


proc initSinglePolynomialTerm(coefficient: Fraction, variable: char, degree: Fraction): PolynomialTerm =
    result.coefficient = coefficient
    result.variable = variable
    result.degree = degree

proc initMultiplePolynomialTerms(variables: openArray[PolynomialTerm],operations: openArray[InfixOperations]): PolynomialTerm =
    assert variables.len div 2 == operations.len, "Too much or too little infix operations provided"
    result.variables = @variables
    result.infixOperations = @operations
    result.degree = 1

proc initPolynomialEquation(variables: openArray[PolynomialTerm], operations: openArray[InfixOperations]): PolynomialEquation =
    assert variables.len div 2 == operations.len, "Too much or too little infix operations provided"
    result.variables = @variables
    result.infixOperations = @operations

converter toPolynomialTerm(a: int): PolynomialTerm =
    result.coefficient = a.toFraction
    result.degree = zeroFraction

converter toPolynomialTerm(a: Fraction): PolynomialTerm =
    result.coefficient = a
    result.degree = zeroFraction

converter toPolynomialEquation(term: PolynomialTerm): PolynomialEquation =
    if term.isMultipleTerms:
        result.variables = term.variables
        result.infixOperations = term.infixOperations
    else:
        result.variables.add(term)

func withVariable(term: PolynomialTerm, variable: char): PolynomialTerm =
    result = term
    result.variable = variable

func inverse(term: PolynomialTerm): PolynomialTerm =
    result = term
    result.degree = term.degree.additiveInverse

    if not term.isMultipleTerms:
        result.coefficient = term.coefficient.multiplicativeInverse

func isSameUnknownTerms(a, b: PolynomialTerm): bool =
    assert a.isMultipleTerms == false and b.isMultipleTerms == false
    result = a.variable == b.variable

func isSameDegree(a, b: PolynomialTerm): bool =
    assert a.isMultipleTerms == false and b.isMultipleTerms == false
    result = a.degree == b.degree

template canAdd(a, b: PolynomialTerm): bool = (isSameUnknownTerms(a, b) and isSameDegree(a, b))
template canMultiply(a, b: PolynomialTerm): bool = isSameUnknownTerms(a, b)

proc `$`(term: PolynomialTerm): string =
    if term.isMultipleTerms:
        result.add(fmt"({$term.variables[0]}) {$term.infixOperations[0]} ({$term.variables[1]})")
        for termIndex in 2..<term.variables.len:
            result.add(fmt"{$term.infixOperations[termIndex+1 div 2]} ({$term.variables[termIndex]})")

        let degreeStr = if term.degree != 1.toFraction: "^" & $term.degree else: ""

        result = fmt"""({result}){degreeStr}"""
    else:
        let
            variable = $(term.variable)

        var
            coefficient, degree: string
            
        if term.coefficient == zeroFraction:
            return ""
        elif term.coefficient == 1.toFraction:
            coefficient = ""
        else:
            coefficient = $term.coefficient

        if term.degree == zeroFraction:
            return coefficient
        elif term.degree == 1.toFraction:
            degree = ""
        else:
            degree = $term.degree
        
        result = coefficient & variable & "^" & degree

func `+`(a, b: PolynomialTerm): PolynomialTerm =
    # assert a.degree == b.degree and a.variable == b.variable, "Unlike terms trying to add"
    if canAdd(a, b):
        result = a 
        result.coefficient = a.coefficient + b.coefficient
    else:
        result.isMultipleTerms = true
        result.variables = @[a, b]
        result.infixOperations = @[Plus]
        result.degree = 1.toFraction

func `-`(a, b: PolynomialTerm): PolynomialTerm =
    # assert a.degree == b.degree and a.variable == b.variable, "Unlike terms trying to subtract"
    if canAdd(a, b):
        result = a 
        result.coefficient = a.coefficient - b.coefficient
    else:
        result.isMultipleTerms = true
        result.variables = @[a, b]
        result.infixOperations = @[Minus]
        result.degree = 1.toFraction

func `*`(a, b: PolynomialTerm): PolynomialTerm =
    # assert a.variable == b.variable, "Unlike Terms trying to multiply"
    if canMultiply(a, b):
        result.variable = a.variable
        result.coefficient = a.coefficient * b.coefficient
        result.degree = a.degree + b.degree
    else:
        result.isMultipleTerms = true
        result.variables = @[a, b]
        result.infixOperations = @[Multiplication]
        result.degree = 1.toFraction


func `/`(a, b: PolynomialTerm): PolynomialTerm =
    # assert a.variable == b.variable, "Unlike Terms trying to divide"\
    if canMultiply(a, b):
        result.variable = a.variable
        result.coefficient = a.coefficient / b.coefficient
        result.degree = a.degree - b.degree
    else:
        result.isMultipleTerms = true
        result.variables = @[a, b]
        result.infixOperations = @[Division]
        result.degree = 1.toFraction

func `*`(a: Fraction, b: PolynomialTerm): PolynomialTerm =
    let constantTerm = a.toPolynomialTerm.withVariable(b.variable)
    result = constantTerm * b

func `+`(a: Fraction, b: PolynomialTerm): PolynomialTerm =
    let constantTerm = a.toPolynomialTerm.withVariable(b.variable)
    result = constantTerm + b

func `-`(a: Fraction, b: PolynomialTerm): PolynomialTerm =
    let constantTerm = a.toPolynomialTerm.withVariable(b.variable)
    result = constantTerm - b

func `-`(a: PolynomialTerm, b: Fraction): PolynomialTerm =
    let constantTerm = b.toPolynomialTerm.withVariable(a.variable)
    result = a - constantTerm

template `+`(a: PolynomialTerm, b: Fraction): PolynomialTerm = b + a
template `*`(a: PolynomialTerm, b: Fraction): PolynomialTerm = b * a
template `/`(a: Fraction, b: PolynomialTerm): PolynomialTerm = a * (b.inverse)
template `/`(a: PolynomialTerm, b: Fraction): PolynomialTerm = a * b.multiplicativeInverse



when isMainModule:
    let x = initSinglePolynomialTerm(5.toFraction, 'x', initFrac(1, 2))
    let y = initSinglePolynomialTerm(5.toFraction, 'y', initFrac(1, 4))
    let z = initSinglePolynomialTerm(5.toFraction, 'z', 2.toFraction)
    echo 5 * (x + y)


#TODO, make sure products of a number and multiple polynomial terms works