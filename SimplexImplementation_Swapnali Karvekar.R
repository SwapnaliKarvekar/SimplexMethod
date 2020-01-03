##Extracting data
##First column has row names 
##after first row, next rows has constraints and the coefficients associated with variables
## Next row has the cost associated with each variables
## Next row has order associated with which we have arranged variables in the columns
## User will input data in such a way that non bases will be entered first and then bases 
##First row has column names 
## Second row has b vector
## Next columns has Matrix A
##n is number of variables in the original question
##m is number of constraints

Data1 <- read.csv(file.choose(), row.names = 1, header= TRUE)
View(Data1)

## In this data file, user will put number of constraints in first column and number of variables in the second column
Data2 <- read.csv(file.choose(), header= TRUE)
View(Data2)

## Extracting m and n
m <- Data2[1,1]
n <- Data2[1,2]

## User will comment if the problem is maximization or minimization
direction = "min"

## Creating Matrix A so that values from data can be extracted
A <- matrix(rep(0,((m+2)*(m+n))), nrow = (m+2) , ncol =  (m+n))
## Extracting Matrix A from the data1
i <-1 
for (i in 1:(m+n)){
  for (j in 1:(m+2)){
    A[j,i] <- Data1[j, (i+1)]
  }
}
A

## Extracting Matrix B from the Matrix A
B = A[1:m, (n+1):(n+m)] 
B
## Extracting Matrix N from the Matrix A
N = A[1:m, 1:n]
N

## Extracting Cost Matrix for bases CB from the Matrix A
CB <- matrix(rep(0,m), nrow = 1 , ncol =  m) 
for (i in 1:m){
  CB[1,i] <- A[(m+1),(n+i)]
}
CB

## Creating Cost Matrix for non bases CN 
CN <- matrix(rep(0,n), nrow = 1 , ncol =  n) 

## Extracting Cost Matrix for non bases CN from the Matrix A
for (i in 1:n){
  CN[1,i] <- A[(m+1),i]
}
CN

## Creating Matrix to keep track of which variables are bases and which are not
## This will be created from Data1 file where user entered the order of variables in the last row
##This matrix will be updated after each iteration
C <- matrix(rep(0,m+n), nrow = 1 , ncol =  m+n)

#### Extracting Order of variables from the Matrix A
i <- 1
for (i in (1:(m+n))){
  C[1,i] <- A[(m+2),i]
}
C
## Extracting Matrix b from the Matrix A
b <-Data1[1:m, 1] 
b
## Creating Matrix Aj(Non bases entering bases) & Ak(Bases leaving the bases) for non bases bases CN
Aj <- matrix(rep(0,m), nrow = m , ncol =  1)
Ak <- matrix(rep(0,m), nrow = m , ncol =  1)
done = FALSE

## Setting interation to 1
p <- 1
while (!done){
  
##Calculating Inverse of bases
  library(matlib)
  BI <- inv (B)
  BI

##Calculating Reduced cost(RC) 
## Doing it step by step as R shows error if done simultaneously
  Z1 <-BI%*%N
  Z2 <- CB%*%Z1
  RC <-CN -Z2

  RC = matrix(RC, nrow =1 )
  RC
  write.csv(RC, "rc.csv")

  
##Creating index t to check if all or some of the reduced cost are postive for minimization problem or negative for maximization problem
## If t= n, the algorithm will be stopped and value of z will be extracted
## If t is not equal to n, another iteration will be performed
  t <- 0
  i <- 1
## Calculating value of t by checking values of RC
  if(direction == "min"){
    for (i in 1:n){
      if (RC [1,i] > 0){
        t <- t +1
      }else{
        i <- i +1
      }
    }
  }else{
    for (i in 1:n){
      if (RC [1,i] < 0){
        t <- t +1
      }else{
        i <- i +1
      }
    }
  }
  t
  
## Determining whether to end the algorith or go for another iteration
## If it is stopped will get optimal value, values of variables(XB), Matrix C(To know which are bases and which are not), number of iterations
## In the else function( if we decide to go for another iteration)
    if (t == n){
    done = TRUE
    X <- (BI%*%b)
    print(X)
    print(C)
    OPtimal=CB%*%X 
    print(OPtimal)
    Iterations <- p
    print(Iterations)
    print(done)
  }else {
    p <- p+1
    i<-1
    j <- (1)
## It will check values of RC and extract Aj(Matrix associated with variable entering the bases) from non bases depending upon minimization of maximization of problem
## e is the order number of variable which will be entering the bases
if(direction == "min"){
      for (j in 1:n){
        if (RC[1, j] < 0){
          for (i in 1:m){
            Aj[i,1] <- N[i,j]
            e <- j
          }
        }
      }
    }else{
      for (j in 1:n){
        if (RC[1, j] > 0){
          for (i in 1:m){
            Aj[i,1] <- N[i,j]
            e <- j
          }
        }
      }
      
    }
    
## Reduced cost of variable selected for entering the bases    
    RC[1, e]
## A matrix associated with variable entering the bases
    Aj
## Order of the matrix entering the bases     
    ## e is calulated so that we can update matrix c and can perform exchanging between bases and non bases    e
    
## Calculating dB 
    dB <- -(BI%*%Aj)
    dB
## Calculating XB(Solution)
    X <- (BI%*%b)
    X
    
## Checking if the problem is unbounded
## Calculating how many number of rows(d) in dB are positive
## If d is equal to m, the problem is unbounded
    d <- 0
    for (i in 1:m){
      if (dB [i,1] >0){
        d <- d+1
      }
    }
    
## If d is equal to m, algorithm will be stopped 
## If not, varible leaving the bases will be calculated 
    if (d == m){
        Solution = print("Unbounded")
        done = TRUE
      }else{
## Calculating which variable is leaving the bases by checking minimum value of lamda(l)
    i <- 1
    for (i in 1:m){
      if (dB [i,1] <0){
        l <- (-(X[i,1]/dB[i]))
      }
    }
    l
    i  

    i<-1
    k<- 0
    for (i in 1:m){
      if (dB[i,1] < 0){
        a <- (-(X[i,1]/dB[i, 1]))
        if (a > 0){
          if ( a < l){
            l <- a
            k <- i
          }
        }
      }
    }
    l
    
  ## k is the order of variable which will be leaving the bases
  ## k is calulated so that we can update matrix c and can perform exchanging between bases and non bases  
    for (i in 1:m){
      if ((-(X[i,1]/dB[i, 1])) == l){
      k <- i
      }
    }
    k
  ## Extracting matrix Ak, which will leave the bases
    for (i in 1:m){
      Ak[i,1] <- B[i, k]
    }

    Ak
  ## Performing exchange between bases and non bases 
  ## First Bases are updated to new values
  ## Then Non bases are updated using matrix Ak
    i <- 1
    for (i in 1:m){
      B[i,k] <- N[i, e]
      N[i,e] <- Ak[i, 1]
    }
    N
    B
    ## Performing exchange between cost matrix associated with bases and non bases 
    ## First value of cost basis associated with variable leaving the bases is assigned to ck
    ## Then Bases are updated to new values
    ## Then Non bases are updated using ck
    ck <- CB[1, k]
    ck
    CB[1, k] <- CN[1, e]
    CN[1, e] <- ck
    CB
    
    ## Updating order matric C to keep the track of variables which are changing from bases to non bases
    ## It is performed in a similar way as updating cost matrix

    f <- C[1, (k)]
    f
    C[1, k] <- C[1, m+e]
    C[1, k]
    C[1, (m+e)] <- f
    C
    print(C)
      }
  }
} 
print(OPtimal)
print(C)
print(Iterations)
View(rc.csv)