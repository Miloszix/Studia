## zad1.1
## zad1.2
function silnia(x)
    a = 1
    
    for i in 1:x
        a *= i
    end
    
    return a
    end
    silnia(5)
   
    ## zad1.3
    function parzystosc(x)
        if x % 2 == 0
            return true
        else
            return false
        end
    end
    
    parzystosc(4)
    parzystosc(5)
    
    ## zad1.4
    function pierwsza(x)
        a = 0
    
        for i in 1:x
            if x % i == 0
            continue
        else
            break
        i += i
        i = a
        end
        
        if a == 0
    
        end
    ##zad 1.7
    function triangle(x)
        a = 1
        for i in 1:x
            a *= 0.75
            i += i
        end
        
        return a
        end
    
        triangle(3)