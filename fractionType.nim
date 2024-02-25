import math

type Fraction = tuple
    numerator: int
    denominator: int

proc initFrac(numerator = 0, denominator = 1):  Fraction =
    result.numerator = numerator
    result.denominator = denominator

converter toFloat(a: Fraction): float = a.numerator / a.denominator

converter toFraction(a: int): Fraction = initFrac(numerator = a)

template zeroFraction(): Fraction = initFrac()

proc multiplicativeInverse(a: Fraction): Fraction =
    if a.numerator == 0:
        result = zeroFraction
    else:
        result = initFrac(a.denominator, a.numerator)

proc additiveInverse(a: Fraction): Fraction =
    if a.numerator == 0:
        result = zeroFraction
    else:
        result = initFrac(a.numerator * -1, a.denominator)

proc `$`(a: Fraction): string =
    if a.denominator == 1: 
        result = $a.numerator #reduce noise
    else:
        result = $a.numerator & "/" & $a.denominator

proc toProperFraction(a: Fraction): Fraction =
    #This doesnt actually make it a "Proper Fraction", it just reduces both numerator and denominator into the lowest values while still retaining the scale. E.g `2/4` is converted to `1/2`
    #Also this cures the fraction of any `div by zero` madness

    let gcd = gcd(a.numerator, a.denominator)

    if (a.denominator < 0):
        let 
            nonNegativeDenominator = a.denominator * -1
            newNumerator = a.numerator * -1
        return initFrac(newNumerator div gcd, nonNegativeDenominator div gcd)
    else:
        if a.denominator == 0:
            return zeroFraction
        else:
            return initFrac(a.numerator div gcd, a.denominator div gcd)


proc `*` (a, b: Fraction): Fraction =
    defer: result = result.toProperFraction

    if (a == zeroFraction) or (b == zeroFraction):
        result = zeroFraction
    else:
        result.numerator = a.numerator * b.numerator
        result.denominator = a.denominator * b.denominator

proc `*`(a: int, b: Fraction): Fraction =
    defer: result = result.toProperFraction

    result = a.toFraction * b

proc `+` (a, b: Fraction): Fraction = 
    defer: result = result.toProperFraction

    if (a == additiveInverse(b)):
        result = zeroFraction
    else:
        result.numerator = (a.numerator*b.denominator) + (b.numerator * a.denominator)
        result.denominator = a.denominator * b.denominator

proc `/`(a, b: Fraction): Fraction =
    defer: result = result.toProperFraction

    result = a * b.multiplicativeInverse

proc `/`(a: int, b: Fraction): Fraction =
    defer: result = result.toProperFraction

    result = a * b.multiplicativeInverse

proc `/`(a: Fraction, b: int):Fraction =
    defer: result = result.toProperFraction

    result = a * b.toFraction.multiplicativeInverse

# template `div` (a, b: Fraction): Fraction = a / b
# template `div` (a, b: int): Fraction = a / b

template `-`(a, b: Fraction): Fraction =
    a + (-1 * b)

template `-`(a: Fraction, b: int): Fraction =
    a + (-1 * b.toFraction)

template `-`(a: int, b: Fraction): Fraction =
    a.toFraction + (-1 * b)

template `+=`(a, b: Fraction) = a = a + b
template `*=`(a, b: Fraction) = a = a * b


proc lcm(a, b: Fraction): Fraction = 
    if a == b:
        return a
    else:
        let product = a * b
        let lcm = lcm(product.numerator, product.denominator)
        return lcm.toFraction