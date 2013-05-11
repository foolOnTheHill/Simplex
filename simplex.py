######################################################################################################
################################################# SIMPLEX ############################################
######################################################################################################
########################Author : George H. A. de Oliveira (ghao@cin.ufpe.br)##########################
######################################################################################################
# -------------------------------------------------------------------------------------------------- #
# The point here is not efficiency, but details and understanding. All the calculations are printed  #
# and this way I wish to make easiest to learn Linear Algebra and its procedures. Have fun.          #
# -------------------------------------------------------------------------------------------------- #
######################################################################################################

import math

def printMatrix(M, decimals = 2) :
    """Imprime uma matriz na forma usual."""
    for x in range(len(M)) :
        s = ""
        for y in range(len(M[x])) :
            if type(M[x][y]) == type(1.0) and M[x][y].is_integer() == False:
                s = s + "\t" + (str(M[x][y])[:str(M[x][y]).index('.') + (decimals + 1)]) #Formatação para duas casas decimais
            else :
                s = s + "\t" + str(M[x][y])
        print(s)

def transpose(M) :
    """Retorna a transposta de M, i.e., M^t ."""
    N = []
    for x in range(len(M[0])) :
        t = []
        for y in range(len(M)) :
            t.append(M[y][x])
        N.append(t)
    return N

def transposeVector(v) :
    """Retorna o vetor v transposto, i.e., v^t ."""
    t = []
    for x in range(len(v[0])) :
        for y in range(len(v)) :
            t.append(v[y][x])
    return t

def sumMatrix(M, N) :
    """Retorna a matriz resultante da soma de M e N."""
    if len(M) == len(N) and len(M[0]) == len(N[0]):
        R = []
        for x in range(len(M)) :
            t = []
            for y in range(len(M[0])) :
                    t.append(M[x][y] + N[x][y])
            R.append(t)
        return R
    else :
        print("As matrizes devem possuir a mesma dimensão!")

def subtractMatrix(M, N) :
    """Retorna a matriz resultante de (M-N)."""
    if len(M) == len(N) and len(M[0]) == len(N[0]) :
        R = []
        for x in range(len(M)) :
            t = []
            for y in range(len(M[0])) :
                t.append(M[x][y] - N[x][y])
            R.append(t)
        return R
    else :
        print("As matrizes devem possuir a mesma dimensão!")

def multiplyMatrix(M, N) :
    """Retorna a matriz resultante de MxN , nessa ordem."""
    if len(M[0]) == len(N) :
                    R = []
                    for i in range(len(M)) : #Iteração pelas linhas de M
                        L = [] #Linhas da Matrix Resultado
                        for j in range(len(N[0])) : #Iteração pelas colunas de N
                            s = 0 #Variável temporária para cada termo 'Rij'
                            for k in range(len(M[i])) : #Iteração pelas colunas das linhas 'i' de M
                                s += (M[i][k]*N[k][j])
                            L.append(s)
                        R.append(L)
                    return R
    else :
                     print("Dimensões Incorretas! O número de colunas de M deve ser igual ao número de linhas de N.")

def multiplyVectorScalar(k, v) :
    """Retorna o vetor resultante da multiplicação de k por v, i.e., k*v ."""
    z = []
    for x in range(len(v)) :
        z.append(k*v[x])
    return z

def divideVectorScalar(k, v) :
    """Retorna o vetor resultante da multiplicação de 1/k por v, i.e., (1/k)*v ."""
    z = []
    for x in range(len(v)) :
        if v[x] == 0 :
            z.append(0.0)
        else :
            z.append(v[x]/float(k))
    return z

def multiplyMatrixScalar(k, M) :
    """Altera M para a resultante de k*M."""
    for x in range(len(M)) :
            M[x] = multiplyVectorScalar(k, M[x])

def sumVectors(u, v) :
    """Retorna o vetor resultante de u+v."""
    if len(u) == len(v) :
        t = []
        for x in range(len(u)) :
            t.append(u[x] + v[x])
        return t
    else :
        print("Os vetores devem possuir a mesma dimensão!")

def subtractVectors(u, v) :
    """Retorna o vetor resultante de u-v."""
    if len(u) == len(v) :
        t = []
        for x in range(len(u)) :
            t.append(u[x] - v[x])
        return t
    else :  
        print("Os vetores devem possuir a mesma dimensão!")

def scalarProduct(u, v) :
    """Retorna o escalar resultante de <u, v> ou u.v ."""
    t = transpose(u)
    r = multiplyMatrix(t, v)[0][0]
    return r

def vectorLength(v) :
    """Retorna a norma ou módulo do vetor v."""
    return math.sqrt(scalarProduct(v, v))

def perm(M) :
    """Realiza as permutações necessárias para o funcionamento correto dos procedimentos de fatoração de matrizes, retornando o número de permutações realizadas."""
    n = 0 #Num de permutações feitas -> usado no cálculo do determinante
    firstLine = False #Verificar se ocorre troca logo na primeira linha

    for x in range(len(M)) :
        """ Procura pelo índice da coluna do pivô """
        for y in range(len(M[x])) :
            if M[x][y] != 0 :
                k = y
                break
            
        """ Procura por um termo não nulo nas próximas linhas
            e numa coluna anterior ao pivô em questão """
        for i in range(x+1, len(M)) :
            for j in range(0, k) :
                if M[i][j] != 0 :
                    if x == 0 : #Troca de primeira linha
                        firstLine = True
                    n += 1
                    t = M[i]
                    M[i] = M[x]
                    M[x] = t
                    print("Trocando a linha " + str(x + 1) + " pela linha " + str(i + 1) + ":")
                    printMatrix(M)
                    print("")
                    break
                
            if firstLine == True:
                break
                
        if firstLine == True :
            break
        
    return n
                
def gaussianElimination(M):
    """Realiza a Eliminação de Gauss sobre a matriz M, alterando-a para a sua forma 'Triangular Superior'."""
    printMatrix(M)
    print("")
    
    n = 0 #Num de permutções -> usado no cálculo do determinante
    for count in range(len(M)) :
        n += perm(M)
        
        """ Busca pelo pivô da linha correspondente """
        pivotColumn = None
        for x in range(len(M[count])) :
            if M[count][x] != 0 :
                pivotColumn = x
                break
        """#########################################"""

        """ Subtração das linhas seguintes """
        for i in range(count + 1, len(M)) :

            if pivotColumn == None : break #Linha nula
            
            if M[count][pivotColumn] != 0 :
                m = M[i][pivotColumn]/float(M[count][pivotColumn])
                if m != 0 :
                    M[i] = sumVectors(multiplyVectorScalar(-m, M[count]), M[i])
                    print("Linha " + str(i + 1) + " menos " + str(m) + " vezes a linha " + str(count + 1) + ": ")
                    printMatrix(M)
                    print("")
        
    return n

def gaussJordanFactorization(M) :
    """Fatora a matriz M para a sua forma escalonada."""
    gaussianElimination(M)
    
    for count in range(len(M) - 1, -1, -1) :
        
        """ Busca pelo pivô da linha correspondente """
        pivotColumn = None
        for x in range(len(M[count])) :
            if M[count][x] != 0 :
                pivotColumn = x
                break
        """#########################################"""
        
        if pivotColumn != None : #Se possuir pivot...
            if M[count][pivotColumn] != 1 :
                print("Dividindo a linha " + str(count + 1) + " por " + str(M[count][pivotColumn]) + ".")
                M[count] = divideVectorScalar(M[count][pivotColumn], M[count])
                printMatrix(M)
                print("")
    
        """ Subtração das linhas seguintes """
        for i in range(count - 1, -1, -1) :

            if pivotColumn == None : break #Linha nula -> pass
                
            if M[count][pivotColumn] != 0 :
                m = M[i][pivotColumn]/float(M[count][pivotColumn])
                if m != 0 :
                    M[i] = sumVectors(multiplyVectorScalar(-m, M[count]), M[i])
                    print("Linha " + str(i + 1) + " menos " + str(m) + " vezes a linha " + str(count + 1) + ": ")
                    printMatrix(M)
                    print("")
        """#################################"""
    return M

def inverse(M) :
    """Retorna a matriz inversa, caso exista, da matriz M dada, i.e., M^-1 ."""
    if len(M) != len(M[0]) :
        print("A matriz deve ser quadrada!")
    else :
        """##################################"""
        W = []
        for x in range(len(M)) :
            t = []
            for y in range(len(M[0])) :
                t.append(M[x][y])
            W.append(t)
        """##################################"""

        singular = False
        terminado = False

        linhas = len(M)
        colunas = len(M[0])
        
        while singular == False  and terminado == False:
            
            printMatrix(W)
            print("")
            print("Adicionando a Matriz Identidade à direita da matriz original: ")
                
            for i in range(len(W)) :
                t = []
                for j in range(len(W[i])) :
                    if i == j :
                        t.append(1)
                    else :
                        t.append(0)
                W[i].extend(t)

            printMatrix(W)
            print("")
                
            for count in range(linhas) :
                perm(W)
                    
                """ Busca pelo pivô da linha correspondente """
                pivotColumn = None
                for x in range(colunas) :
                    if W[count][x] != 0 :
                        pivotColumn = x
                        break
                """#########################################"""

                """ Subtração das linhas seguintes """
                for i in range(count + 1, linhas) :

                    if pivotColumn == None :
                        singular = True
                        break #Linha nula
                        
                    if W[count][pivotColumn] != 0 :
                        m = W[i][pivotColumn]/float(W[count][pivotColumn])
                        if m != 0 :
                            W[i] = sumVectors(multiplyVectorScalar(-m, W[count]), W[i])
                            print("Linha " + str(i + 1) + " menos " + str(m) + " vezes a linha " + str(count + 1) + ": ")
                            printMatrix(W)
                            print("")
                                
            for count in range(linhas - 1, -1, -1) :
                    
                """ Busca pelo pivô da linha correspondente """
                pivotColumn = None
                for x in range(colunas) :
                    if W[count][x] != 0 :
                        pivotColumn = x
                        break
                """#########################################"""
                    
                if pivotColumn != None : #Se possuir pivot...
                    if W[count][pivotColumn] != 1 :
                        print("Dividindo a linha " + str(count + 1) + " por " + str(W[count][pivotColumn]) + ".")
                        W[count] = divideVectorScalar(W[count][pivotColumn], W[count])
                        printMatrix(W)
                        print("")
                    
                """ Subtração das linhas seguintes """
                for i in range(count - 1, -1, -1) :
                        
                    if pivotColumn == None :
                        singular = True
                        break #Linha nula -> pass
                            
                    if W[count][pivotColumn] != 0 :
                        m = W[i][pivotColumn]/float(W[count][pivotColumn])
                        if m != 0 :
                            W[i] = sumVectors(multiplyVectorScalar(-m, W[count]), W[i])
                            print("Linha " + str(i + 1) + " menos " + str(m) + " vezes a linha " + str(count + 1) + ": ")
                            printMatrix(W)
                            print("")

                terminado = True              
                "#################################"""
                    
        if singular :
            print("Ops... A matriz é singular .'. não possui inversa.")
        else :
            I = []
            for k in range(len(M)) :
                index = int(len(W[k])/2)
                I.append(W[k][index:])
            print("A inversa é : ")
            printMatrix(I)
            print("")
            return I

def dependence(V) :
    """Testa a Independência Linear do conjunto de vetores V dado."""

    """ Cria uma cópia do conjunto original para nao alterá-lo """
    W = []
    singular = False
    terminado = False
    
    for x in range(len(V)) :
        W.append(transposeVector(V[x]))
    """########################################################"""

    r = "O conjunto " + str(W) + " é " 
    
    while singular == False and terminado == False:
        
        for count in range(len(W)) :
            perm(W)
            
            """ Busca pelo pivô da linha correspondente """
            pivotColumn = None
            for x in range(len(W[count])) :
                if W[count][x] != 0 :
                    pivotColumn = x
                    break
            """#########################################"""

            """ Subtração das linhas seguintes """
            for i in range(count + 1, len(W)) :

                if pivotColumn == None :
                    singular = True
                    break
                
                if W[count][pivotColumn] != 0 :
                    m = float(W[i][pivotColumn]/W[count][pivotColumn])
                    if m != 0 :
                        W[i] = sumVectors(multiplyVectorScalar(-m, W[count]), W[i])
                        
        for count in range(len(W) - 1, -1, -1) :
            
            """ Busca pelo pivô da linha correspondente """
            pivotColumn = None
            for x in range(len(W[count])) :
                if W[count][x] != 0 :
                    pivotColumn = x
                    break
            """#########################################"""
            
            if pivotColumn != None : #Se possuir pivot...
                if W[count][pivotColumn] != 1 :
                    W[count] = divideVectorScalar(W[count][pivotColumn], W[count])
                    
            """ Subtração das linhas seguintes """
            for i in range(count - 1, -1, -1) :

                if pivotColumn == None :
                    singular = True
                    break
                
                if W[count][pivotColumn] != 0 :
                    m = float(W[i][pivotColumn])/W[count][pivotColumn]
                    if m != 0 :
                        W[i] = sumVectors(multiplyVectorScalar(-m, W[count]), W[i])
                        
        terminado = True
        """#################################"""

    if singular == True :
        r += "Linearmente Dependente(LD)."
    else :
        r += "Linearmente Independente(LI)."

    print(r)
    return not(singular)

def determinant(M) :
    """Calcula o determinante da matriz quadrada M."""
    if len(M) != len(M[0]) :
        print("A matriz deve ser quadrada!")
    else :
        n = gaussianElimination(M)
        d = 1
        s = "O determinante é igual ao produtório: \ndet = [(-1)^" + str(n) + "] * ";
        for x in range(len(M)) : #Determinante == 'Produto dos pivots'
            if x == len(M) - 1 :
                s += str(M[x][x]) + "."
            else :
                s += str(M[x][x]) + " * "
            d *= M[x][x]
        d = ((-1)**n) * d
        s += "\ndet = " + str(d)
        print(s)
        return d

def projectionMatrix(M) :
    """Retorna a matriz de projeção no Espaço Coluna da matriz M."""
    
    print("A Matriz de Projeção no Espaço Coluna de M é dada por : P = M*[(M^t*M)^(-1)]*M^t")
    print("M : ")
    printMatrix(M)
    #------------------------------#
    print("\nM^t : ")
    end = transpose(M)
    printMatrix(end)
    #------------------------------#
    print("\nA^t*A : ")
    middle = multiplyMatrix(end, M)
    printMatrix(middle)
    #------------------------------#
    print("\nCálculo da Inversa : ")
    middle = inverse(middle)
    print("(M^t*M)^(-1) : ")
    printMatrix(middle)
    #------------------------------#
    if middle != None :
        print("\nP : ")
        n = multiplyMatrix(middle, end)
        P = multiplyMatrix(M, n)
        printMatrix(P)
        return P

def gramSchmidt(V, normalizar) :
    """Transforma a base V numa Base Ortonormal."""
    ##############################
    def inner(u, v) :
        r = 0
        for x in range(len(u)) :
            r += (u[x] * v[x])
        return r
    ##############################
    def subtract(u, v) :
        r = []
        for x in range(len(u)) :
            r.append(u[x] - v[x])
        return r
    ##############################
    Q = []
    for i in range(len(V)) :
        q = []
        for k in range(len(V[i])) :
            q.append(V[i][k])
        #------------------------#
        for j in range(i) :
            s = multiplyVectorScalar((inner(V[i], Q[j])/float(inner(Q[j], Q[j]))), Q[j])
            q = subtract(q, s)
        #------------------------#
        Q.append(q)
    ##############################
    if normalizar == True :
        for x in range(len(Q)) :
            l = math.sqrt(inner(Q[x], Q[x]))
            Q[x] = divideVectorScalar(l, Q[x])
    ###############################
    return Q
    
def volume(P) :
    """Retorna o volume da figura formada pelo conjunto P de pontos."""
    N = []
    for x in range(len(P)) :
        t = []
        for y in range(len(P[x])) :
            t.append(P[x][y])
        t.append(1)
        N.append(t)
    dimension = len(P[0])
    d = determinant(N)
    if d != None :
        v = (1.0/dimension)*d
        if v < 0 : v *= -1
        if dimension == 2 :
            print("A área é : A = " + (str(v)[:4]) + " u.a .")
        else :
            print("O volume é : V = " + (str(v)[:4]) + " u.v .")
        return v

#Examples(Just press F5!) :#
M = [[1, 2, 1, 5, 1], [2, 4, 1, 7, 8], [1, 1, 0, 0, 1]]
N = [[1, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 0], [0, 0, 1]]
P = multiplyMatrix(M, N)
print("M :")
printMatrix(M)
print("N :")
printMatrix(N)
print("P = MxN:")
printMatrix(P)
print("")
gaussJordanFactorization(M)
determinant(P)
