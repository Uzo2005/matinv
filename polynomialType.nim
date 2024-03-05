import ./fractionType

type 
    PolynomialTerm* = tuple
        coefficient: Fraction
        variable: char
        degree: Fraction
        # derivative: Fraction

    PolynomialEquation* = seq[PolynomialTerm] #all polynomial equations implicitly equal 0

proc `$`(term: PolynomialTerm): string =
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
        return "1"
    elif term.degree == 1.toFraction:
        degree = ""
    else:
        degree = $term.degree
    
    result = coefficient & variable & "^" & degree

proc initPolynomialTerm(coefficient: Fraction, variable: char, degree: Fraction): PolynomialTerm =
    discard


when isMainModule:
    discard